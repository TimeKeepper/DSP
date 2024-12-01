clear; % Clear all variables from the workspace
clc; % Clear the command window

N = 1024; % Define the length of the signal
xn = 0:N-1; % Create a vector from 0 to N-1

% Generate two cosine signals
x1 = cos(pi*xn/4); % Cosine signal with frequency pi/4
x2 = cos(pi*xn/4) + cos(pi*xn/8); % Sum of two cosine signals with frequencies pi/4 and pi/8

% Repeat the first 8 samples of x1 to fill the length N
x1 = repmat(x1(1:8), 1, N/8);
% Repeat the first 16 samples of x2 to fill the length N
x2 = repmat(x2(1:16), 1, N/16);

% Perform FFT and plot the results for x1
fft_and_plot(x1, xn, [8, 16]);
figure; % Create a new figure window
% Perform FFT and plot the results for x2
fft_and_plot(x2, xn, [8, 16]);

% Function to perform FFT and plot the results
function fft_and_plot(x, xn, Ns)
    for i = 1:length(Ns) % Loop over the different FFT lengths
        N = Ns(i); % Get the current FFT length
        Xk = fft(x, N); % Compute the FFT of the signal with length N
        subplot(2, 1, i); % Create a subplot for the current FFT length
        stem(xn(1:N)/N, abs(Xk)); % Plot the magnitude of the FFT
        xlim([0 1]); % Set the x-axis limits from 0 to 1
    end
end
