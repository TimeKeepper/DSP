% Clear workspace and command window
clear
clc

% Define the index variable xn from 0 to 1023
xn = 0:1023;

% Initialize three signals x1, x2, and x3 with zeros
x1 = zeros(1,1024);
x2 = zeros(1,1024);
x3 = zeros(1,1024);

% Define the signal x1: set to 1 for indices 0 to 3
x1((xn<=3)) = 1;

% Define the signal x2:
% - Set to xn + 1 for indices 0 to 3
% - Set to 8 - xn for indices 4 to 7
x2((xn<=3)) = xn((xn<=3)) + 1;
x2((xn>=4) & (xn<=7)) = 8 - xn((xn>=4) & (xn<=7));

% Define the signal x3:
% - Set to 4 - xn for indices 0 to 3
% - Set to xn - 3 for indices 4 to 7
x3((xn<=3)) = 4 - xn((xn<=3));
x3((xn>=4) & (xn<=7)) = xn((xn>=4) & (xn<=7)) - 3;

% Store the signals and their titles in cell arrays
signals = {x1, x2, x3};
titles = {'x1', 'x2', 'x3'};

% Define an array of FFT lengths to be used
N_values = [1024, 8, 16];

% Loop through each signal
for i = 1:3
    % Create a new figure for each signal
    figure
    % Loop through each FFT length
    for j = 1:3
        % Get the current FFT length
        N = N_values(j);
        % Compute the FFT of the current signal with the current length
        Xk = fft(signals{i}, N);
        % Create a subplot for the current FFT result
        subplot(3, 1, j);
        % Plot the magnitude of the FFT result
        if N == 1024
            % Use plot for N=1024
            plot(xn(1:N)/N, abs(Xk));
        else
            % Use stem for N=8 and N=16
            stem(xn(1:N), abs(Xk));
        end
        % Title the subplot with the corresponding signal name and FFT length
        title(['FFT of ', titles{i}, ' with N=', num2str(N)]);
    end
end
