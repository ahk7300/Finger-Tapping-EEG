% Load the EEG data as a matrix
filename = 'BCI_10AT_10IT.csv'; % Replace with your file name
data = readmatrix(filename); % Read the data

% Define columns for time and selected EEG channels
time_column = 1; % Assuming time is in the first column
channel_columns = 2:9; % Adjust these indices based on your EEG channels

% Extract time and EEG signals
time = data(:, time_column); % Time in seconds
eeg = data(:, channel_columns); % EEG signals for selected channels

% Sampling frequency and parameters
fs = 250; % Sampling frequency in Hz (adjust if different)
low_cutoff = 0.5; % Low-frequency cutoff in Hz
high_cutoff = 50; % High-frequency cutoff in Hz

% Design a bandpass filter
[b, a] = butter(4, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');
filtered_eeg = filtfilt(b, a, eeg); % Apply zero-phase filtering

% Define baseline and task segments (adjust based on your experiment)
baseline_start = 0; % Baseline start time in seconds
baseline_end = 2;   % Baseline end time in seconds
task_start = 3;     % Task start time in seconds
task_end = 5;       % Task end time in seconds

% Convert times to indices
baseline_idx = (time >= baseline_start & time <= baseline_end);
task_idx = (time >= task_start & time <= task_end);

% Band definitions
alpha_band = [8 12]; % Alpha band
beta_band = [13 30]; % Beta band

% Filter for alpha and beta bands
[b_alpha, a_alpha] = butter(4, alpha_band / (fs / 2), 'bandpass');
alpha_filtered = filtfilt(b_alpha, a_alpha, filtered_eeg);

[b_beta, a_beta] = butter(4, beta_band / (fs / 2), 'bandpass');
beta_filtered = filtfilt(b_beta, a_beta, filtered_eeg);

% Calculate power during baseline and task for alpha and beta bands
baseline_power_alpha = mean(alpha_filtered(baseline_idx, :).^2, 1);
task_power_alpha = mean(alpha_filtered(task_idx, :).^2, 1);

baseline_power_beta = mean(beta_filtered(baseline_idx, :).^2, 1);
task_power_beta = mean(beta_filtered(task_idx, :).^2, 1);

% Calculate ERD/ERS percentage for alpha and beta bands
erd_ers_alpha = ((task_power_alpha - baseline_power_alpha) ./ baseline_power_alpha) * 100;
erd_ers_beta = ((task_power_beta - baseline_power_beta) ./ baseline_power_beta) * 100;

% Display ERD/ERS results
disp('ERD/ERS for Alpha Band (%):');
for i = 1:size(eeg, 2)
    fprintf('Channel %d: %.2f%%\n', i, erd_ers_alpha(i));
end

disp('ERD/ERS for Beta Band (%):');
for i = 1:size(eeg, 2)
    fprintf('Channel %d: %.2f%%\n', i, erd_ers_beta(i));
end

% Plot time-domain signals for alpha and beta bands
figure;
for i = 1:size(eeg, 2)
    subplot(size(eeg, 2), 1, i);
    plot(time, alpha_filtered(:, i), 'b'); hold on;
    plot(time, beta_filtered(:, i), 'r');
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(['Channel ', num2str(i), ': Alpha (Blue) and Beta (Red) Signals']);
    legend('Alpha Band', 'Beta Band');
    grid on;
end
saveas(gcf, 'time_domain_alpha_beta_signals.png');

% Frequency-domain analysis (FFT) for alpha and beta bands
figure;
for i = 1:size(eeg, 2)
    n = length(alpha_filtered(:, i));
    f = (0:n-1) * (fs / n); % Frequency vector

    fft_alpha = abs(fft(alpha_filtered(:, i)) / n);
    fft_beta = abs(fft(beta_filtered(:, i)) / n);

    subplot(size(eeg, 2), 1, i);
    plot(f(1:n/2), fft_alpha(1:n/2), 'b'); hold on;
    plot(f(1:n/2), fft_beta(1:n/2), 'r');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title(['Channel ', num2str(i), ': Alpha (Blue) and Beta (Red) Frequency Domain']);
    legend('Alpha Band', 'Beta Band');
    grid on;
end
saveas(gcf, 'frequency_domain_alpha_beta_signals.png');

% PSD for alpha and beta bands using Welch's method
figure;
for i = 1:size(eeg, 2)
    [psd_alpha, f_psd_alpha] = pwelch(alpha_filtered(:, i), hamming(256), 128, 256, fs);
    [psd_beta, f_psd_beta] = pwelch(beta_filtered(:, i), hamming(256), 128, 256, fs);

    subplot(size(eeg, 2), 1, i);
    plot(f_psd_alpha, 10*log10(psd_alpha), 'b'); hold on;
    plot(f_psd_beta, 10*log10(psd_beta), 'r');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title(['Channel ', num2str(i), ': Alpha (Blue) and Beta (Red) PSD']);
    legend('Alpha Band', 'Beta Band');
    grid on;
end
saveas(gcf, 'psd_alpha_beta.png');

disp('Processing complete. ERD/ERS values and plots saved.');
