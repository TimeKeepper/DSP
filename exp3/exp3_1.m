clc; % 清除命令窗口
clear; % 清除工作区变量

% 定义输入信号 x
% x = [0, 1, 2, 3, 4, 5, 6, 7];
x = [2, 3, 4, 5, 6];
N = length(x); % 获取输入信号的长度

p = 0:N-1; % 定义时间索引
m = 0:N-1; % 定义频率索引

% 预分配频域信号 X 的大小
X = zeros(1, N);
% 计算离散傅里叶变换 (DFT)
for k = 1:N
    for n = 1:N
        % 计算 X(k) 的值
        X(k) = X(k) + x(n) * exp(-1i * 2 * pi / N) .^ ((n-1) * (k-1));
    end
end

% 预分配逆变换信号 x1 的大小
x1 = zeros(1, N);
% 计算逆离散傅里叶变换 (IDFT)
for n = 1:N
    for k = 1:N
        % 计算 x1(n) 的值
        x1(n) = x1(n) + X(k) * exp(1i * 2 * pi / N) .^ ((n-1) * (k-1)) / N;
    end
end

% 绘制原始信号 x(n)
subplot(2,2,1);
stem(p, x);
title('(a) x(n)');

% 绘制逆变换结果 x1(n) 的幅度
subplot(2,2,2);
stem(p, abs(x1));
title('(b) IDFT 结果 x1(n)');

% 绘制频域信号 X(k) 的幅度谱
subplot(2,2,3);
stem(m, abs(X));
title('(c) X(k)的幅度谱');

% 绘制频域信号 X(k) 的相位谱
subplot(2,2,4);
stem(m, angle(X));
title('(d) X(k)的相位谱');
