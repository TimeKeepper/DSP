clear
clc

% Parameters
f1 = 50;
f2 = 100;
fs = 1000;
Ts = 1/fs;
N = 160;
k = 0.5;
t = 0:Ts:(N-1)*Ts;

% Generate signals
xt = cos(2*pi*f1*t) + cos(2*pi*f2*t) + k*randn(1, length(t));
xn = xt(1:N);

% Plot time domain waveform
subplot(3,1,1)
plot(t, xt);
title("时域波形");

% Compute and plot power spectral density with noise
Xk = fft(xn);
P = abs(Xk).^2 / N;
subplot(3,1,2)
stem((0:N-1)/N, P);
title("功率谱密度");

% Generate noise-free signal and compute power spectral density
xt_no_noise = cos(2*pi*f1*t) + cos(2*pi*f2*t);
xn_no_noise = xt_no_noise(1:N);
Xk_no_noise = fft(xn_no_noise);
P_no_noise = abs(Xk_no_noise).^2 / N;
subplot(3,1,3)
stem((0:N-1)/N, P_no_noise);
title("无噪声下功率谱密度");