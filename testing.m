% Load EEG data
filename = 'BCI_10AT_10IT.csv'; % Replace with the actual path if needed
data = readmatrix(filename); % Load data matrix

% Define time and EEG channels
time_column = 1; % Assuming time is in the first column
channel_columns = 2:9; % Adjust based on EEG channels

time = data(:, time_column); % Extract time
eeg = data(:, channel_columns); % Extract EEG signals

% Sampling frequency and parameters
fs = 250; % Sampling frequency in Hz
low_cutoff = 0.5; % Low cutoff frequency in Hz
high_cutoff = 50; % High cutoff frequency in Hz

% Bandpass filter design
[b, a] = butter(4, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');
filtered_eeg = filtfilt(b, a, eeg);

% Define intervals
baseline_interval = [0, 2]; % Baseline time interval
actual_tap_intervals = [0, 10; 15, 25]; % Actual tap intervals
imagined_tap_intervals = [10, 15; 25, 30]; % Imagined tap intervals

% Convert intervals to indices
baseline_idx = time >= baseline_interval(1) & time <= baseline_interval(2);

% Corrected loop-based indexing for actual and imagined taps
actual_tap_idx = false(size(time)); % Initialize as false
for j = 1:size(actual_tap_intervals, 1)
    actual_tap_idx = actual_tap_idx | (time >= actual_tap_intervals(j, 1) & time <= actual_tap_intervals(j, 2));
end

imagined_tap_idx = false(size(time)); % Initialize as false
for j = 1:size(imagined_tap_intervals, 1)
    imagined_tap_idx = imagined_tap_idx | (time >= imagined_tap_intervals(j, 1) & time <= imagined_tap_intervals(j, 2));
end

% Band definitions
alpha_band = [8 12]; % Alpha band
beta_band = [13 30]; % Beta band

% Filter for alpha and beta bands
[b_alpha, a_alpha] = butter(4, alpha_band / (fs / 2), 'bandpass');
alpha_filtered = filtfilt(b_alpha, a_alpha, filtered_eeg);

[b_beta, a_beta] = butter(4, beta_band / (fs / 2), 'bandpass');
beta_filtered = filtfilt(b_beta, a_beta, filtered_eeg);

% Calculate power during baseline, actual tap, and imagined tap for alpha and beta bands
baseline_power_alpha = mean(alpha_filtered(baseline_idx, :).^2, 1);
actual_tap_power_alpha = mean(alpha_filtered(actual_tap_idx, :).^2, 1);
imagined_tap_power_alpha = mean(alpha_filtered(imagined_tap_idx, :).^2, 1);

baseline_power_beta = mean(beta_filtered(baseline_idx, :).^2, 1);
actual_tap_power_beta = mean(beta_filtered(actual_tap_idx, :).^2, 1);
imagined_tap_power_beta = mean(beta_filtered(imagined_tap_idx, :).^2, 1);

% Calculate ERD/ERS percentages for alpha and beta bands
erd_ers_alpha_actual = ((actual_tap_power_alpha - baseline_power_alpha) ./ baseline_power_alpha) * 100;
erd_ers_alpha_imagined = ((imagined_tap_power_alpha - baseline_power_alpha) ./ baseline_power_alpha) * 100;

erd_ers_beta_actual = ((actual_tap_power_beta - baseline_power_beta) ./ baseline_power_beta) * 100;
erd_ers_beta_imagined = ((imagined_tap_power_beta - baseline_power_beta) ./ baseline_power_beta) * 100;

% Display ERD/ERS results
disp('ERD/ERS for Alpha Band (%):');
for i = 1:size(eeg, 2)
    fprintf('Channel %d - Actual Tap: %.2f%%, Imagined Tap: %.2f%%\n', ...
            i, erd_ers_alpha_actual(i), erd_ers_alpha_imagined(i));
end

disp('ERD/ERS for Beta Band (%):');
for i = 1:size(eeg, 2)
    fprintf('Channel %d - Actual Tap: %.2f%%, Imagined Tap: %.2f%%\n', ...
            i, erd_ers_beta_actual(i), erd_ers_beta_imagined(i));
end

% Save ERD/ERS results to a file
results_table = table((1:size(eeg, 2))', erd_ers_alpha_actual', erd_ers_alpha_imagined', ...
                      erd_ers_beta_actual', erd_ers_beta_imagined', ...
                      'VariableNames', {'Channel', 'Alpha_Actual', 'Alpha_Imagined', ...
                                        'Beta_Actual', 'Beta_Imagined'});
writetable(results_table, 'ERD_ERS_results.csv');
disp('ERD/ERS results saved to "ERD_ERS_results.csv".');

% Plot filtered EEG signals with tap markers
figure;
for i = 1:size(eeg, 2)
    subplot(size(eeg, 2), 1, i);
    plot(time, filtered_eeg(:, i), 'b'); hold on;

    % Add markers for actual taps
    for j = 1:size(actual_tap_intervals, 1)
        x = actual_tap_intervals(j, :);
        fill([x(1) x(2) x(2) x(1)], ...
             [min(filtered_eeg(:, i)) min(filtered_eeg(:, i)) max(filtered_eeg(:, i)) max(filtered_eeg(:, i))], ...
             'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Green for actual taps
    end

    % Add markers for imagined taps
    for j = 1:size(imagined_tap_intervals, 1)
        x = imagined_tap_intervals(j, :);
        fill([x(1) x(2) x(2) x(1)], ...
             [min(filtered_eeg(:, i)) min(filtered_eeg(:, i)) max(filtered_eeg(:, i)) max(filtered_eeg(:, i))], ...
             'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Red for imagined taps
    end

    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(['Channel ', num2str(i), ': Filtered EEG with Tapping Markers']);
    grid on;
end
saveas(gcf, 'filtered_eeg_tapping_markers.png');
disp('Plot with tapping markers saved as "filtered_eeg_tapping_markers.png".');
