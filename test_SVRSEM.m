clc;
clear;

% add path
run('\package\irt/setup');  
run('\package\KER_v0.11/setup'); 

N=256;               % 图像大小
N2=N^2;
imgsiz = [N N]; 
disp('Generating phantom image');
load('slice_x_90_256×256_pet_lesion.mat')
I = model;
Image=reshape(I,N2,1);        % 产生实时图像数据

theta=linspace(0,360,211);
theta=theta(1:210);   % 投影角度

%% ========= 产生投影数据 ========= %%
disp('产生投影数据');
P_num=249;           % 探测器通道个数
M=P_num*length(theta);                  % 投影射线的总条数

%% ========= 获取系统矩阵 ========= %%
disp('获取系统矩阵')
sysflag = 1;  % 1: using Fessler's IRT; 
              % 0: using your own system matrix G 
disp('--- Generating system matrix G ...')
if sysflag
    % require Fessler's IRT matlab toolbox to generate a system matrix
    ig = image_geom('nx', imgsiz(1), 'ny', imgsiz(2), 'fov', 33.3);

    % field of view
    ig.mask = ig.circ(16, 16) > 0;

    % system matrix G
    prjsiz = [249 210];
    sg = sino_geom('par', 'nb', prjsiz(1), 'na', prjsiz(2), 'dr', 70/prjsiz(1), 'strip_width', 'dr');
    G  = Gtomo2_strip(sg, ig, 'single', 1);
    Gopt.mtype  = 'fessler';
    Gopt.ig     = ig;
    Gopt.imgsiz = imgsiz;
else
    Gfold = ''; % where you store your system matrix G;
    load([Gfold, 'G']);	% load matrix G
    Gopt.mtype  = 'matlab';
    Gopt.imgsiz = imgsiz;
end


  x0 = I;
  x0 = reshape(x0, [], 1);
  for m = 1:size(x0,2)
      proj(:,m) = proj_forw(G, Gopt, x0(:,m));
  end

  % attenuation factor
  u = zeros(1,N2);
  u(model ~= 0) = 0.0098;
  ai = exp(-repmat(proj_forw(G, Gopt, u(:)), 1));  

  % background (randoms and scatters)
  ri = repmat(mean(ai.*proj,1)*0.2,[size(proj,1) 1]);  % 20% uniform background

  % total noiseless projection
  y0 = ai.*proj + ri;

  % count level
  count = 250e6; 

  % normalized sinograms
  cs = count / sum(y0(:));
  y0 = cs * y0;
  ri = cs * ri;
  yi=poissrnd(y0);
  ni = ai*cs; % multiplicative factor
  yi=reshape(yi,M,1);

  % ----- 提取系统矩阵 ----- %
  load("A.mat")

  % ----- 提取权重矩阵 ----- %
  load('weightMatrix_256.mat')
  w = weightMatrix;


  % % ------ 加载最优(BSREM) ------ % %
  load("BSREM_model_250_10.mat")
  x_star = BSREM_model_250_10;

%% ========= 划分子集 ========= %%
disp('划分子集')
L=10;             % 子集个数
T=length(theta)/L;        % 每个子集包含的角度数
theta_seq=reshape(1:length(theta),L,T);       % 产生角度编号矩阵
W_index=zeros(L,T*P_num);            % 包含L个子集的射线的行号
for i=1:L
    temp=P_num*theta_seq(i,:);
    for j=1:T
        W_index(i,(j-1)*P_num+1 :j*P_num)=temp(j)-P_num+1:temp(j);
    end
end

load("indexes_1.mat")
load('indexes_2.mat')

%% = = = = = = = 利用BSREM_SVRG算法进行重建 = = = = = = = %%
disp('开始进行BSREM_SVRG重建')
F0=ones(N2,1);             % 初始图像向量
irt_numS=10;               % 迭代次数

err5=zeros(irt_numS+1,1);
[F5,err5,t5,Phi5,mse5,mse5_lesion1,mse5_lesion2,suboptimality_k] = medfuncSVRSEM(W_index,N,F0,yi,ni,ri,x0,L,G,A,w,Gopt,irt_numS,err5,Image,indexes_1,indexes_2,x_star);
F5=reshape(F5,N,N);

%% ========= 仿真结果显示 ========= %%
% 绘制CPU运行时间与误差之间的图像
figure(1);
plot(0:0.1:irt_numS*(2*L)/10,err5,'b-','LineWidth',1);
xlabel('epochs');
ylabel('NRMSD');

figure(2);
plot(0:0.1:(irt_numS*(2*L)/10)-0.1,suboptimality_k,'b-','LineWidth',1);

figure(3);
imshow(F5),xlabel('Proposed','FontSize',15);
colormap('jet')
clim([0 1]);
title(sprintf('MSE = %5.2fdB',mse5),FontSize=13);




