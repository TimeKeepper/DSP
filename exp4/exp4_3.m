% Clear the workspace and command window
clear
clc

% Sampling frequency
fs = 64;

% Define the range of n
n = -1024:1024;

% Generate the signal xn as a sum of three cosine waves
xn = cos(8*pi*n/fs) + cos(16*pi*n/fs) + cos(20*pi*n/fs);

% Define the different values of N for the FFT
N_values = [16, 32, 64];

% Loop over each value of N
for i = 1:length(N_values)
    N = N_values(i); % Current FFT length
    Xk = fft(xn, N); % Compute the N-point FFT of the signal
    Xk = fftshift(Xk); % Shift zero frequency component to the center
    F = fs/N; % Frequency resolution
    subplot(3,1,i); % Create a subplot for each N
    stem((-N/2:N/2-1)*F, abs(Xk)); % Plot the magnitude of the FFT
    title(['N = ', num2str(N)]); % Title indicating the value of N
end
