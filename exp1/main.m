% 基本矩阵运算
M = [1 2 3; 4 5 6; 7 8 9]
M_trans = M'

M = [1+2j 3+4j]
M'
M.'

A = [1 2; 3 4]
A + 3
A = [1 1 1; 1 1 1]
B = [3 6; 2 4; 5 10]
A*B
A*6

A = [1 2; 3 4]
B = [5 6; 7 8]
A / B
A./ B
A = [2 7 3; 6 5 4; 1 8 9]
inv(A)

% 基本信号生成

% 单位冲激序列
n = 1 : 20
n0 = 8
data = [zeros(1, n0-1) 1 zeros(1, length(n) - n0)]
stem(n, data, 'k')
title('单位冲激序列')
xlabel('n')
ylabel('x(n)')

% 矩阵序列
n1 = 15
un = [zeros(1, n0 - 1) ones(1, n1 - n0 + 1) zeros(1, length(n) - n1)]
stem(n, un, 'k')
title('单位冲激序列')
xlabel('n')
ylabel('x(n)')

% 实指数序列
n = 0: 50
a = 0.9
xn = a.^n
stem(n, xn, 'k')
title('实指数序列')
xlabel('n')
ylabel('x(n)')

% 正弦序列及其采样
n = 0: 100
a = 1.5
f = 0.02
x_sin = a * sin(2*pi*f*n)
subplot(2, 1, 1)
plot(n, x_sin)
title('正弦序列')
xlabel('n')
ylabel('x(n)')
grid on
un = 1 : 33 : length(n) % 此处修改频率
xn = x_sin([un]);
subplot(2, 1, 2)
stem(un - 1, xn)
title('采样正弦序列')
xlabel('n')
ylabel('x(n)')
grid on

% 随机序列
x1 = rand(1, 5)
figure(1)
stem(x1)
a = 4;
b = 6;
x2 = a + (b - a) * rand(1, 5)
figure(2)
stem(x2)

% 序列的移位
x = [1 4  3 2];
nx = 1:length(x);
m = 3;
ny = nx + m;
y = x;
subplot(211);
stem(nx, x);
axis([min(min(nx), min(ny))-1, max(max(nx), max(ny))+1, min(x)-1, max(x)+1]);
xlabel('n')
ylabel('x(n)')
title('原始序列')
grid on
subplot(212);
stem(ny, y);
axis([min(min(nx), min(ny))-1, max(max(nx), max(ny))+1, min(x)-1, max(x)+1]);
xlabel('n')
ylabel('y(n)')
title('移位后的序列')
grid on

% 序列的求和与乘积
clc;
clear;
x1 = [1, 3, 5, 7, 6, 4, 2, 1];
ns1 = -3;
x2 = [4, 0, 2, 1, -1, 3];
ns2 = 1;
nf1 = ns1 + length(x1) - 1;
nf2 = ns2 + length(x2) - 1;
n1 = ns1 : nf1;
n2 = ns2 : nf2;
n = min(ns1, ns2) : max(nf1, nf2);
y1 = zeros(1, length(n));
y2 = y1;
y1(find((n >= ns1) & (n <= nf1) == 1)) = x1;
y2(find((n >= ns2) & (n <= nf2) == 1)) = x2;
ya = y1 + y2;
ym = y1.* y2;

subplot(221);
stem(n1, x1, '.');
ylabel('x1(n)');
grid;

subplot(223);
stem(n2, x2, '.')
xlabel('n');
ylabel('x2(n)');
grid;

subplot(222);
stem(n, ya, '.')
ylabel('y1(n) + y2(n)');
grid;

subplot(224);
stem(n, ym, '.')
xlabel('n');
ylabel('y1(n) * y2(n)');
grid;

% 求序列的能量（功率）
clc;
clear;
xn = [1 + 1i, 2 + 2i, 3 + 3i, 4 + 4i];
E = sum(abs(xn).^2)

n = 0:100;
a = 1.5;
f = 0.02;
x_sin = a * sin(2*pi*f*n);
E_sin = sum(abs(x_sin).^2)
P = E_sin / n(end)

% 序列的卷积运算
clc;
clear;
x1 = [1, 3, 5, 7];
ns1 = -3;
x2 = [4, 0, 2, 1];
ns2 = 1;
nf1 = ns1 + length(x1) - 1;
nf2 = ns2 + length(x2) - 1;
nx1 = ns1 : nf1;
nx2 = ns2 : nf2;
ny1 = ns1 + ns2;
ny2 = nf1 + nf2;
ny = ny1 : ny2;
y = conv(x1, x2);

figure(1);
subplot(121);
stem(nx1, x1);
xlabel('n');
ylabel('x1(n)');
subplot(122);
stem(nx2, x2);
xlabel('n');
ylabel('x2(n)');
grid;

figure(2);
stem(ny, y)
xlabel('n');
ylabel('y(n)');
grid;

% % 数字信号处理 byd实验不一样这个不用做
% clc;
% clear;
% [yy, fs] = audioread('final.wav');
% y = yy(:, 1);
% sound(y, fs);
% py = fft(y);
% figure(1);
% subplot(1, 2, 1), plot(1:length(y), y)
% subplot(1, 2, 2), plot(1:length(py)/2, py(1:length(py)/2))
% title('原始音频')

% [b, a] = butter(9, 0.25, 'high');
% yh = filter(b, a, y);
% sound(yh, fs)
% pyh = fft(yh);
% figure(2);
% subplot(1, 2, 1), plot(1:length(yh), yh)
% subplot(1, 2, 2), plot(1:length(pyh)/2, pyh(1:length(pyh)/2))
% title('高通滤波后')
% audiowrite('final_high.wav', yh, fs)

% [b, a] = butter(9, 0.05, 'low');
% yl = filter(b, a, y);
% sound(yl, fs)
% pyl = fft(yl);
% figure(3);
% subplot(1, 2, 1), plot(1:length(yl), yl)
% subplot(1, 2, 2), plot(1:length(pyl)/2, pyl(1:length(pyl)/2))
% title('低通滤波后')
% audiowrite('final_low.wav', yl, fs)
