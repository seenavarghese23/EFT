% Battery parameters
V_nominal = 3.7;  % Nominal voltage in V
Cn = 4.5 * 3600;  % Nominal capacity in As (Ampere-seconds)
R0 = 0.002;      % Ohmic resistance
R1 = 0.003;      % First RC pair resistance
C1 = 1000e-3;    % First RC pair capacitance
R2 = 0.00522;    % Second RC pair resistance
C2 = 1000e-3;    % Second RC pair capacitance

% Simulation parameters
t_end = 3600;   % Simulation time in seconds
dt = 1;         % Time step in seconds
t = 0:dt:t_end; % Time vector
N = length(t);  % Number of simulation steps

% Load SOC time data
addpath 'C:\Users\DELL\Documents\Master Thesis EFT\github\Old as 17-07-2024\Test'
load SOC_time_data1234.mat;
SOC_true = data1;

% Generate simulated current profile (e.g., constant discharge)
I_true = -5 * ones(1, N);  % 5A discharge current

% Generate true voltage measurements based on the battery model
V_true = ocv(SOC_true) - I_true * R0 - I_true * R1 .* (1 - exp(-t ./ (R1 * C1))) - I_true * R2 .* (1 - exp(-t ./ (R2 * C2)));

% Define noise levels for comparison
voltage_noise_levels = [0.01, 0.05, 0.1];  % Standard deviations for voltage noise
current_noise_levels = [0.001, 0.05, 0.1];  % Standard deviations for current noise

% Initialize arrays to store results for different noise types
SOC_est_results_aekf = zeros(length(voltage_noise_levels), length(current_noise_levels), N);

% Helper function to generate pink noise
function x = pinknoise(N)
    white_noise = randn(1, N);
    x = filter([0.021 0.044 0.021], [1 -2.347 1.977 -0.626], white_noise);
end

% Run simulations for different noise levels and types
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        voltage_noise_std = voltage_noise_levels(v);
        current_noise_std = current_noise_levels(c);
        
        % Generate noisy measurements
        V_meas = V_true + voltage_noise_std * randn(1, N);
        I_meas = I_true + current_noise_std * randn(1, N);

        % Initialize AEKF
        x_aekf = [100; 2; 2];  % Initial state [SOC; V_RC1; V_RC2]
        P = diag([0.01^1, 0.01^1, 0.01^1]);  % Initial state covariance
        Q = diag([10e-7, 10e-5, 10e-5]);     % Process noise covariance
        R = diag([max(voltage_noise_std^2, 1e-10), max(current_noise_std^2, 1e-10)]);  % Measurement noise covariance

        % AEKF-specific parameters for adaptive noise tuning
        lambda = 1;  % Forgetting factor
        alpha = 0.1; % Adaptation rate for measurement noise covariance
        gamma = 0.01; % Innovation scaling factor

        % AEKF loop
        for k = 1:N
            % Prediction step
            [x_pred, F] = state_transition(x_aekf, I_meas(k), dt, Cn, R1, C1, R2, C2);
            P_pred = F * P * F' + Q;
            
            % Update step
            z = [V_meas(k); I_meas(k)];
            [z_pred, H] = measurement_model(x_pred, I_meas(k), R0, R1, R2);
            y = z - z_pred;
            S = H * P_pred * H' + R;
            K = P_pred * H' / S;
            x_aekf = x_pred + K * y;
            P = (eye(3) - K * H) * P_pred;

            % Innovation-based adaptive tuning of R
            R = lambda * R + (1 - lambda) * (y * y') / (H * P_pred * H' + gamma);

            % Store results
            SOC_est_results_aekf(v, c, k) = x_aekf(1);
        end
    end
end

% Save results after all simulations
save('SOC_est_results_aekf.mat', 'SOC_est_results_aekf', 'voltage_noise_levels', 'current_noise_levels');

% Extract and save SOC vs. time data for AEKF
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        soc_vs_time_aekf = [t'/3600, squeeze(SOC_est_results_aekf(v, c, :))];
        filename_aekf = sprintf('SOC_vs_time_aekf_Vnoise_%.3f_Inoise_%.3f.csv', voltage_noise_levels(v), current_noise_levels(c));
        writematrix(soc_vs_time_aekf, filename_aekf);
    end
end

% Plot results for SOC estimation using AEKF
figure;
plot_idx = 1;
legend_labels = {};  % Initialize empty cell array for legend labels
colors = lines(length(voltage_noise_levels) * length(current_noise_levels));  % Generate distinct colors dynamically

for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t/3600, squeeze(SOC_est_results_aekf(v, c, :)), 'Color', colors(plot_idx, :), 'LineWidth', 1.5);
        hold on;
        legend_labels{plot_idx} = sprintf('V noise: %.3f, I noise: %.3f', voltage_noise_levels(v), current_noise_levels(c));  % Ensure string format
        plot_idx = plot_idx + 1;
    end
end

plot(t/3600, SOC_true, 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'True SOC';

xlabel('Time (hours)');
ylabel('SOC');
legend(legend_labels, 'Location', 'southwest');
title('SOC Estimation with AEKF');
grid on;

% Helper functions
function [x_next, F] = state_transition(x, i, dt, Cn, R1, C1, R2, C2)
    % State transition model
    x_next = zeros(3, 1);
    x_next(1) = x(1) - (i * dt) / Cn;  % SOC update
    x_next(2) = exp(-dt / (R1 * C1)) * x(2) + R1 * (1 - exp(-dt / (R1 * C1))) * i;  % V_RC1 update
    x_next(3) = exp(-dt / (R2 * C2)) * x(3) + R2 * (1 - exp(-dt / (R2 * C2))) * i;  % V_RC2 update
    
    % State transition Jacobian matrix
    F = [1, 0, 0;
         0, exp(-dt / (R1 * C1)), 0;
         0, 0, exp(-dt / (R2 * C2))];
end

function [z, H] = measurement_model(x, i, R0, R1, R2)
    % Measurement model
    z = [ocv(x(1)) - x(2) - x(3) - i * R0; i];
    
    % Measurement Jacobian matrix
    H = [docv_dsoc(x(1)), -1, -1;
         0, 0, 0];
end

function v = ocv(soc)
    % Open circuit voltage as a function of SOC
    v = 3.7 + 0.3 * soc;
end

function dv = docv_dsoc(soc)
    % Derivative of OCV with respect to SOC
    dv = 0.3;
end
figure;

% White noise subplot
subplot(3, 1, 1);
plot_idx = 1;
legend_labels = {};  % Reset legend labels

for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t/3600, squeeze(SOC_est_results_white(v, c, :)), 'Color', colors(plot_idx, :), 'LineWidth', 1.5);
        hold on;
        legend_labels{plot_idx} = sprintf('V noise: %.3f, I noise: %.3f', voltage_noise_levels(v), current_noise_levels(c));
        plot_idx = plot_idx + 1;
    end
end

plot(t/3600, SOC_true, 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'True SOC';

xlabel('Time (hours)');
ylabel('SOC');
legend(legend_labels, 'Location', 'southwest');
title('SOC Estimation with White Noise');
grid on;

% Gaussian noise subplot
subplot(3, 1, 2);
plot_idx = 1;
legend_labels = {};  % Reset legend labels

for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t/3600, squeeze(SOC_est_results_gaussian(v, c, :)), 'Color', colors(plot_idx, :), 'LineWidth', 1.5);
        hold on;
        legend_labels{plot_idx} = sprintf('V noise: %.3f, I noise: %.3f', voltage_noise_levels(v), current_noise_levels(c));
        plot_idx = plot_idx + 1;
    end
end

plot(t/3600, SOC_true, 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'True SOC';

xlabel('Time (hours)');
ylabel('SOC');
legend(legend_labels, 'Location', 'southwest');
title('SOC Estimation with Gaussian Noise');
grid on;

% Pink noise subplot
subplot(3, 1, 3);
plot_idx = 1;
legend_labels = {};  % Reset legend labels

for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t/3600, squeeze(SOC_est_results_pink(v, c, :)), 'Color', colors(plot_idx, :), 'LineWidth', 1.5);
        hold on;
        legend_labels{plot_idx} = sprintf('V noise: %.3f, I noise: %.3f', voltage_noise_levels(v), current_noise_levels(c));
        plot_idx = plot_idx + 1;
    end
end

plot(t/3600, SOC_true, 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'True SOC';

xlabel('Time (hours)');
ylabel('SOC');
legend(legend_labels, 'Location', 'southwest');
title('SOC Estimation with Pink Noise');
grid on;
% Calculate deviations
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        % White Noise
        abs_diff_white = abs(squeeze(SOC_est_results_white(v, c, :)) - SOC_true);
        %deviation_white(v, c) = mean(abs_diff_white);  % Use mean for average deviation
        
        % Alternatively, for RMSD (uncomment to use RMSD instead)
        deviation_white(v, c) = sqrt(mean((abs_diff_white).^2));

        % Gaussian Noise
        abs_diff_gaussian = abs(squeeze(SOC_est_results_gaussian(v, c, :)) - SOC_true);
        deviation_gaussian(v, c) = sqrt(mean((abs_diff_gaussian).^2));
        
        % Pink Noise
        abs_diff_pink = abs(squeeze(SOC_est_results_pink(v, c, :)) - SOC_true);
        deviation_pink(v, c) = sqrt(mean((abs_diff_pink).^2));
    end
end

% Create table for deviations
row_names = {};  % Initialize empty cell array for row names
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        row_names{end+1} = sprintf('V_noise = %.3f, I_noise = %.3f', voltage_noise_levels(v), current_noise_levels(c));
    end
end

% Combine deviation results into one table
T = table(deviation_white(:), deviation_gaussian(:), deviation_pink(:), ...
    'VariableNames', {'Deviation_White', 'Deviation_Gaussian', 'Deviation_Pink'}, ...
    'RowNames', row_names);

% Display the table
disp(T);

% Save the table as a CSV file
writetable(T, 'SOC_deviation_table.csv', 'WriteRowNames', true);
