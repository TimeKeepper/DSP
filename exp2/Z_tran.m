% % 有限长序列的Z变换
% clc;clear
% n = 1 : 20
% n0 = 8
% n1 = 15
% % x = [1 zeros(1, length(n) - 1)] % 单位冲激序列
% % x = [zeros(1, n0 - 1) ones(1, n1 - n0 + 1) zeros(1, length(n) - n1)] % 矩阵序列
% % n = 0: 100
% % a = 1.5
% % f = 0.02
% % x = a * sin(2*pi*f*n) % 正弦序列
% n = 0: 50
% a = 0.9
% x = a.^n % 实指数序列
% zx=-1:0.01:1; zy=zx;
% for n=1:length(zx)
%     for m=1:length(zy)
%         z(n,m)=zx(n)+i*zy(m);
%         Zt_tmp=0;
%         for ii=1:length(x)
%             Zt_tmp=Zt_tmp+x(ii)*((z(n,m))^(-(ii - 1)));
%         end
%         Zt(n,m)=Zt_tmp;
%     end
% end
% mesh(zx,zy,abs(Zt))
% title("实指数序列的Z变换")
% hold on;  % 保持当前图形

% 无限长序列的Z变换

clc;clear
zx=-1:0.01:1; zy=zx;
for n=1:length(zx)
    for m=1:length(zy)
        z(n,m)=zx(n)+i*zy(m);
        Zt_tmp=1/(1-0.9*(z(n, m)^-1)); % x[n] = 2^n
        Zt(n,m)=Zt_tmp;
    end
end
figure(1);
mesh(zx,zy,abs(Zt))
title("无限长实指数序列的Z变换")
hold on;  % 保持当前图形

% 定义圆的参数
theta = linspace(0, 2*pi, 100);  % 角度
r = 1;  % 圆的半径
z_circle = 10;  % 圆在 Z 轴的高度

% 计算圆的坐标
x_circle = r * cos(theta);  % 圆的 X 坐标
y_circle = r * sin(theta);  % 圆的 Y 坐标

% 绘制圆
plot3(x_circle, y_circle, z_circle * ones(size(x_circle)), 'r-', 'LineWidth', 2);
