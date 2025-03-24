clear;clc;

load("BSREM_model_250_10.mat")
load("SVRSEM_model_250_10.mat")
load("MLEM_model_250_10.mat")
load("OSEM_model_250_10.mat")
load("SDPBSREM_model_250_10.mat")




% 计算两个图像之间的绝对误差
% error_map_BSREM = abs(BSREM_model_250_10 - SVRSEM_model_250_10);
% error_map_MLEM = abs(BSREM_model_250_10 - MLEM_model_250_10);
error_map_MLEM = BSREM_model_250_10 - MLEM_model_250_10;
error_map_OSEM = BSREM_model_250_10 - OSEM_model_250_10;
error_map_SVRSEM = BSREM_model_250_10 - SVRSEM_model_250_10;
error_map_SDPBSREM = BSREM_model_250_10 - SDPBSREM_model_250_10;


% 显示误差图
figure(1);
hold on;

imshow(error_map_MLEM,[]);  % 绘制误差图
colormap('jet');     % 设置颜色映射为热图（红色表示较大的误差）
colorbar;            % 显示颜色条
% title('MLEM');  % 图表标题
% xlabel('Column Index');  % X轴标签
% ylabel('Row Index');     % Y轴标签
zoom(3); 

hold off;

figure(2);
hold on;

imshow(error_map_OSEM,[]);  % 绘制误差图
colormap('jet');     % 设置颜色映射为热图（红色表示较大的误差）
colorbar;            % 显示颜色条
% title('OSEM');  % 图表标题
% xlabel('Column Index');  % X轴标签
% ylabel('Row Index');     % Y轴标签
zoom(3); 

hold off;


figure(3);
hold on;

imshow(error_map_SDPBSREM,[]);  % 绘制误差图
colormap('jet');     % 设置颜色映射为热图（红色表示较大的误差）
colorbar;            % 显示颜色条
% title('SDP-BSREM');  % 图表标题
% xlabel('Column Index');  % X轴标签
% ylabel('Row Index');     % Y轴标签
zoom(3); 

hold off;


figure(4);
hold on;

imshow(error_map_SVRSEM,[]);  % 绘制误差图
colormap('jet');     % 设置颜色映射为热图（红色表示较大的误差）
colorbar;            % 显示颜色条
% title('Proposed');  % 图表标题
% xlabel('Column Index');  % X轴标签
% ylabel('Row Index');     % Y轴标签
zoom(3); 

hold off;

