% Load the EEG data as a matrix
filename = 'baseline.csv';
data = readmatrix(filename); % Use readmatrix for numeric-only files

% Define columns for time and selected channels (manually specify indices)
time_column = 1; % Assuming time is in the first column
channel_columns = [2, 3, 4, 5, 6, 7, 8, 9]; % Adjust these indices for your channels

% Extract time and EEG signals
time = data(:, time_column); % First column as time
eeg = data(:, channel_columns); % Selected columns for EEG

% Parameters for preprocessing
fs = 250; % Sampling frequency in Hz (adjust based on your data)
low_cutoff = 0.5; % Low-frequency cutoff in Hz
high_cutoff = 50; % High-frequency cutoff in Hz

% Bandpass filter design
[b, a] = butter(4, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');
filtered_eeg = filtfilt(b, a, eeg);

% Time-domain plot for each channel
figure;
for i = 1:size(eeg, 2)
    subplot(size(eeg, 2), 1, i);
    plot(time, filtered_eeg(:, i));
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(['Filtered EEG Signal: Channel ', num2str(i)]);
    grid on;
end
saveas(gcf, 'time_domain_filtered_signals.png');

% Frequency-domain analysis (FFT) for each channel
figure;
for i = 1:size(eeg, 2)
    % Frequency-domain analysis
    n = length(filtered_eeg(:, i));
    f = (0:n-1) * (fs / n); % Frequency vector
    fft_eeg = abs(fft(filtered_eeg(:, i)) / n);

    % Plot FFT result
    subplot(size(eeg, 2), 1, i);
    plot(f(1:n/2), fft_eeg(1:n/2)); % Plot only positive frequencies
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title(['Frequency Domain: Channel ', num2str(i)]);
    grid on;
end
% saveas(gcf, 'frequency_domain_filtered_signals.png');
% 
% disp('Processing complete. Plots saved as PNG files.');
