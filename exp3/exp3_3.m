clc;
clear;

% 定义输入信号 x
x = [0, 1, 2, 3, 4, 5, 6, 7];
N = length(x); % 获取输入信号的长度

p = 0:N-1; % 定义时间索引
w = linspace(-2*pi, 2*pi, 500); % 定义频率范围，从 -2π 到 2π，分成 500 个点

% 初始化频域信号 X
X = zeros(1, 500);

% 计算频域信号 X
for i = 1:500
    for n = 1:N
        X(i) = X(i) + x(n) * exp(-1i * (n-1) * w(i));
    end
end

% 绘制原始信号 x(n)
subplot(3,1,1);
stem(p, x);
title('(a)');
xlabel('n');
ylabel('序列 x(n)');

% 绘制幅度谱
subplot(3,1,2);
plot(w, abs(X));
title('(b)');
xlabel('w');
ylabel('幅度谱');

% 绘制相位谱
subplot(3,1,3);
plot(w, angle(X));
title('(c)');
xlabel('w');
ylabel('相位谱');