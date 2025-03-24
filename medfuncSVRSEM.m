function [F5,err5,t5,Phi5,mse5,mse5_lesion1,mse5_lesion2,suboptimality_k] = medfuncSVRSEM(W_index,N,F0,yi,ni,ri,x0,L,G,A,w,Gopt,irt_numS,err5,Image,indexes_1,indexes_2,x_star)
         
%SVRSEM
% ---------------------------------------------
% 输入参数：
% F0:初始图像向量
% P:投影数据，列向量
% irt_num:迭代次数
% L:子集个数
% ---------------------------------------------
% 输出参数：
% F:输出图像向量
% =============================================%

m = 2*L;

totalElapsedTime = 0; % 用于累积总时间


% set Gopt
Gopt = setGopt(ni, G, Gopt);

% initialization
F5    = max(mean(F0(:))*1e-9,F0(:)); F5(~Gopt.mask) = 0;
epps=1e-8;

err5 = zeros(m,irt_numS);
err5_ini = sqrt(sum((F5-Image).^2))/sqrt(sum((Image).^2));

Phi5 = zeros(m,irt_numS);

MSE = zeros(m,irt_numS);

suboptimality_k = zeros(m,irt_numS);

% 病变区域
lesion1_Image = Image(indexes_1);
lesion2_Image = Image(indexes_2);

t5 = zeros(m,irt_numS);

%% ---------- 计算上界U的值 ---------- %%
for kk=1:L
    pj(:,kk) = A'*ones(size(A,1),1,'double');
    [n,~] = size(A);
    rowMin = nan(n, 1);
    [I,~,S] = find(A);
    [I,K] = sort(I);
    S = S(K);
    markers = [find([1; diff(I)]); numel(I)+1];
    iRows = I(markers(1:end-1));
    for i = 1:numel(iRows)
        s = S(markers(i):(markers(i+1)-1));
        rowMin(iRows(i)) = min(s);
    end
    rowMin(isnan(rowMin)) = 0;
    % Amin(index{ll}) = rowMin;
    Amin = rowMin + epps;
end
D = sum(pj,2);

U = max(double(yi)./Amin);
pj3 = D/L;
for n = 1:size(pj3,2)
    temp = Gopt.ig.embed(pj3(:,n));
    pj_3(:,n) = temp(:);
end
dU = zeros(size(F5,1),1);


t(1,1) = 1;
lamda0 = 0.004;
a = 1/25;

% 在循环前初始化数组
F5_old=F5;


for k=1:irt_numS          % 外循环迭代

    % 在每次迭代开始时启动计时器
    fprintf('迭代 %d / %d\n', k, irt_numS);
    
    F5_tilde = F5_old;

    zb_all = compute_fullgrad(F5_tilde,ni,yi,ri,G,Gopt,W_index,L);
    zb_avg = zb_all/L;
    F5_ol = F5_old;

    for i=1:m

        tic;

        % 随机选择一个子集
            kk = randi(L);

        fprintf('子集 %d / %d\n', kk, L);

        % 更新lamda
        lamda = lamda0/(a*k+1);

        % Nesterov 动量步长
        % 更新 t
        t(k, i+1) = (1 + sqrt(1 + 4 * t(k, i)^2))/2;

        % 更新 alpha
        alpha(k, i) =( 1 + (t(k, i) - 1)/t(k, i+1));

        % 总步长
        p(k,i) = alpha(k,i)*lamda;
   
        % 获取子集需要的系统矩阵和投影数据
        G1 = G(W_index(kk, :), :);
        ni1 = ni(W_index(kk, :), :);
        yi1 = yi(W_index(kk,:),:);
        ri1 = ri(W_index(kk, :), :);

        % % ------- 计算预处理器dU -------%
        pp = F5_ol<U/2;
        dU(pp) = F5_ol(pp)./(pj_3(pp)+epps);
        dU(~pp) = (U-F5_ol(~pp))./(pj_3(~pp)+epps);

        %  ------- 计算数值梯度 -------- %%
        F = reshape(F5_ol,N,N);
        [Fx,Fy]=gradient(F);
        grad_f = sqrt((Fx).^2+(Fy).^2);
        mean_f = mean(F(:));
        grad_f = grad_f(:);
        Mu = max(0.01,grad_f./mean_f);

        V1 = 1.9 ; V2 = 2.05;
        P_v = min(V2,max(mean(Mu)./Mu,V1));

        V_ki = P_v;

        % 随子集变化的预处理器 %
        S_f = V_ki.*dU;
        
        % ------ 计算当前迭代梯度 ------- %
        yb = ni1.*proj_forw(G1, Gopt, F5_ol) +ri1;
        yy = yi1./(yb+epps);
        % yy = yy-1;
        yy(yb==0&yi1==0) = 1;
        zb = proj_back(G1, Gopt, ni1.*(yy-1));

        % ------ 计算之前迭代梯度 ------- %
        yb_tilde = ni1.*proj_forw(G1, Gopt, F5_tilde) +ri1;
        yy_tilde = yi1./(yb_tilde+epps);
        % yy = yy-1;0.096
        yy_tilde(yb_tilde==0&yi1==0) = 1;
        zb_old = proj_back(G1, Gopt, ni1.*(yy_tilde-1));
      
        % 拉普拉斯项的梯度
        F5_diff = F5_ol-F5_ol';
        g_L = sum(w*F5_diff, 2);

        g_svr = (zb - zb_old) + zb_avg;
        % 
        g = g_svr + 0.000001*g_L;


        % % 分块计算法
        % blockSize = 1024;  % 定义块的大小，根据内存限制进行调整
        % n = 65536;  % 矩阵的维度
        % 
        % % 初始化结果向量
        % g = zeros(n, 1);
        % 
        % % 逐块计算
        % for ib = 1:blockSize:n
        %     % 计算当前块的范围
        %     iEnd = min(ib + blockSize - 1, n);
        % 
        %     % 提取 w 的子矩阵
        %     w_block = w(ib:iEnd, :);
        % 
        %     % 计算 F5_ol 的子向量差
        %     F5_diff_block = F5_ol(ib:iEnd) - F5_ol';
        % 
        %     % 逐块进行乘积并累加结果
        %     for j = 1:blockSize:n
        %         jEnd = min(j + blockSize - 1, n);
        % 
        %         % 提取 F5_diff 的子块
        %         F5_diff_subblock = F5_diff_block(:, j:jEnd);
        % 
        %         % 计算并累加结果
        %         g(ib:iEnd) = g(ib:iEnd) + sum(w_block(:, j:jEnd) .* F5_diff_subblock, 2);
        %     end
        % end
        % 
        % g_all = g_svr + 0.000001*g;

        F5_new = F5_ol + p(k,i)*S_f.*g ;


        elapsedTime = toc;
        t5(i,k) = elapsedTime;
        
        err5(i,k) = (sqrt(sum((F5_new-Image).^2))/sqrt(sum((Image).^2)))/err5_ini*100;


        % % 计算次优 % %
        x_k = reshape(F5_new,256,256);
        suboptimality_k(i,k) = norm(x_k - x_star, 'fro');
        

        % 迭代后的病变区域
        lesion1 = F5_new(indexes_1);
        lesion2 = F5_new(indexes_2);
        
        mse5_lesion1(i,k) = 10*log10(sum((lesion1-lesion1_Image).^2)/sum(lesion1_Image.^2));
        mse5_lesion2(i,k) = 10*log10(sum((lesion2-lesion2_Image).^2)/sum(lesion2_Image.^2));

        MSE(i,k) = 10*log10(sum((F5_new-Image).^2)/sum(Image.^2));

        iy = yb>0;
        Phi5(i,k) = sum(yi1(iy).*log(yb(iy))-yb(iy));

        F5_ol = F5_new;

    end

    % 步长相关更新
    t(k+1,1) = t(k,m+1);
    
    F5_old = F5_new;
       
    fprintf('迭代 %d 花费: %.2f 秒\n', k, elapsedTime);
    
    totalElapsedTime = totalElapsedTime + elapsedTime; % 累积总时间

end

suboptimality_k = suboptimality_k(:);
err5 = err5(:); 
err5 = [100;err5];
mse5_lesion1 = mse5_lesion1(:);
mse5_lesion2 = mse5_lesion2(:);
Phi5 = Phi5(:);
Phi5 = [0;Phi5];
t5 = t5(:);
t5 = [0;t5];
fprintf('MBSREM总时间: %.2f 秒\n', totalElapsedTime); % 显示总时间
mse5 = 10*log10(sum((F5_new-Image).^2)/sum(Image.^2)); 

F5 = F5_new;

