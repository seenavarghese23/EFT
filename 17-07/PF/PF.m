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
I_true = 5 * ones(1, N);  % 5A discharge current

% Generate true voltage measurements based on the battery model
V_true = ocv(SOC_true) - I_true * R0 - I_true * R1 .* (1 - exp(-t ./ (R1 * C1))) - I_true * R2 .* (1 - exp(-t ./ (R2 * C2)));

% Define noise levels for comparison
voltage_noise_levels = [0.01, 0.05, 0.1];  % Standard deviations for voltage noise
current_noise_levels = [0.001, 0.05, 0.1];  % Standard deviations for current noise

% Initialize arrays to store results for different noise types
SOC_est_results_white = zeros(length(voltage_noise_levels), length(current_noise_levels), N);
SOC_est_results_gaussian = zeros(length(voltage_noise_levels), length(current_noise_levels), N);
SOC_est_results_pink = zeros(length(voltage_noise_levels), length(current_noise_levels), N);

% Initialize arrays to store average deviations for each noise combination
deviation_white = zeros(length(voltage_noise_levels), length(current_noise_levels));
deviation_gaussian = zeros(length(voltage_noise_levels), length(current_noise_levels));
deviation_pink = zeros(length(voltage_noise_levels), length(current_noise_levels));

% Helper function to generate pink noise
function x = pinknoise(N)
    white_noise = randn(1, N);
    x = filter([0.021 0.044 0.021], [1 -2.347 1.977 -0.626], white_noise);
end

% Number of particles for the particle filter
M = 2000;

% Run simulations for different noise levels and types
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        voltage_noise_std = voltage_noise_levels(v);
        current_noise_std = current_noise_levels(c);
        
        % Generate noisy measurements (White Noise)
        V_meas_white = V_true + voltage_noise_std * randn(1, N);
        I_meas_white = I_true + current_noise_std * randn(1, N);

        % Generate noisy measurements (Gaussian Noise using randn)
        V_meas_gaussian = V_true + voltage_noise_std * randn(1, N);
        I_meas_gaussian = I_true + current_noise_std * randn(1, N);

        % Generate noisy measurements (Pink Noise)
        V_meas_pink = V_true + voltage_noise_std * pinknoise(N);
        I_meas_pink = I_true + current_noise_std * pinknoise(N);

% Initialize particles (SOC, V_RC1, V_RC2) and weights
x_particles = [100 + 5 * randn(1, M); 1 * ones(1, M); 1 * ones(1, M)];  % SOC, V_RC1, V_RC2
weights = ones(1, M) / M;  % Uniform initial weights

% Particle Filter loop for White Noise
for k = 1:N
    % Propagate each particle using the state transition function
    for m = 1:M
        x_particles(:, m) = state_transition(x_particles(:, m), I_meas_white(k), dt, Cn, R1, C1, R2, C2);
    end
    
    % Measurement update and weight calculation
    for m = 1:M
        % Predict the voltage measurement for the current particle
        z_pred = measurement_model(x_particles(:, m), I_meas_white(k), R0, R1, R2);
        
        % Actual noisy measurement (voltage)
        z = V_meas_white(k);
        % Gaussian likelihood function
        weights(m) = weights(m) * exp(-0.5 * (z - z_pred)^2 / voltage_noise_std^2);
    end
    
    % Normalize weights
            weights = weights / sum(weights);
            % Resample particles
            indices = systematic_resample(weights, M);
            x_particles = x_particles(:, indices);
            weights = ones(1, M) / M;
            % Estimate SOC as mean of particles
            SOC_est_results_white(v, c, k) = mean(x_particles(1, :));
    
end


        % Calculate absolute deviation for White Noise
        abs_diff_white = abs(squeeze(SOC_est_results_white(v, c, :)) - SOC_true);
        deviation_white(v, c) = mean(abs_diff_white);  % Store mean deviation

        % PF loop for Gaussian Noise
        x_particles = [100 + 5 * randn(1, M); 2 * ones(1, M); 2 * ones(1, M)];  % SOC, V_RC1, V_RC2
        weights = ones(1, M) / M;  % Uniform initial weights
        for k = 1:N
            % Propagate each particle using state transition
            for m = 1:M
                x_particles(:, m) = state_transition(x_particles(:, m), I_meas_gaussian(k), dt, Cn, R1, C1, R2, C2);
            end
            % Measurement update and weight calculation
            for m = 1:M
                z_pred = measurement_model(x_particles(:, m), I_meas_gaussian(k), R0, R1, R2);
                z = V_meas_gaussian(k);
                weights(m) = weights(m) * exp(-0.5 * sum((z - z_pred).^2) / voltage_noise_std^2);
            end
            % Normalize weights
            weights = weights / sum(weights);
            % Resample particles
            indices = systematic_resample(weights, M);
            x_particles = x_particles(:, indices);
            weights = ones(1, M) / M;
            % Estimate SOC as mean of particles
            SOC_est_results_gaussian(v, c, k) = mean(x_particles(1, :));
        end

        % Particle Filter for Pink Noise
        x_particles = [100 + 5 * randn(1, M); 2 * ones(1, M); 2 * ones(1, M)];
        weights = ones(1, M) / M;
        for k = 1:N
            % Propagate each particle using state transition
            for m = 1:M
                x_particles(:, m) = state_transition(x_particles(:, m), I_meas_pink(k), dt, Cn, R1, C1, R2, C2);
            end
            % Measurement update and weight calculation
            for m = 1:M
                z_pred = measurement_model(x_particles(:, m), I_meas_pink(k), R0, R1, R2);
                z = V_meas_pink(k);
                weights(m) = weights(m) * exp(-0.5 * sum((z - z_pred).^2) / voltage_noise_std^2);
            end
            % Normalize weights
            weights = weights / sum(weights);
            % Resample particles
            indices = systematic_resample(weights, M);
            x_particles = x_particles(:, indices);
            weights = ones(1, M) / M;
            % Estimate SOC as mean of particles
            SOC_est_results_pink(v, c, k) = mean(x_particles(1, :));
        end

        % Calculate absolute deviation for Gaussian and Pink Noise
        abs_diff_gaussian = abs(squeeze(SOC_est_results_gaussian(v, c, :)) - SOC_true);
        deviation_gaussian(v, c) = mean(abs_diff_gaussian);
        
        abs_diff_pink = abs(squeeze(SOC_est_results_pink(v, c, :)) - SOC_true);
        deviation_pink(v, c) = mean(abs_diff_pink);
    end
end

% Plot results or compare deviations for white, gaussian, and pink noise SOC estimates
function ocv_val = ocv(SOC)
    % Open-circuit voltage (OCV) function
    ocv_val = 3.6 + 0.4 * SOC;  % Example linear relationship
end
% Define the state transition function
function x_next = state_transition(x_curr, I, dt, Cn, R1, C1, R2, C2)
    SOC = x_curr(1);   % Current SOC
    V_RC1 = x_curr(2);  % Voltage across the first RC pair
    V_RC2 = x_curr(3);  % Voltage across the second RC pair
    % Update the SOC based on the current
    SOC_next = SOC - (I * dt / Cn);  % SOC decreases as current is drawn
    SOC_next = max(SOC_next, 0);     % Ensure SOC does not drop below 0%
   %disp(SOC_next);
    % Update the RC circuit voltages using first-order RC models
    V_RC1_next = V_RC1 * exp(-dt / (R1 * C1)) + I * R1 * (1 - exp(-dt / (R1 * C1)));
    V_RC2_next = V_RC2 * exp(-dt / (R2 * C2)) + I * R2 * (1 - exp(-dt / (R2 * C2)));
    % Return the updated state vector
    x_next = [SOC_next; V_RC1_next; V_RC2_next];
end
% Define the measurement model
function z = measurement_model(x_curr, I, R0, R1, R2)
    SOC = x_curr(1);
    V_RC1 = x_curr(2);
    V_RC2 = x_curr(3);
    % Open-circuit voltage (OCV) as a function of SOC
    V_OCV = ocv(SOC);  % You can define the OCV function based on SOC lookup or empirical data
    % Predicted voltage based on current and RC components
    V_pred = V_OCV - I * R0 - V_RC1 - V_RC2;
    % Return voltage as the measurement
    z = V_pred;
end
% Systematic Resampling function
function indices = systematic_resample(weights, M)
    % Normalize weights
    weights = weights / sum(weights);
    % Generate cumulative sum of weights
    cumsum_weights = cumsum(weights);
    % Generate positions spaced evenly in [0,1] with random start offset
    positions = (0:(M-1)) / M + rand / M;
    % Initialize variables
    indices = zeros(1, M);
    i = 1;  % Pointer for positions
    j = 1;  % Pointer for cumsum_weights
    % Loop over each position and resample
    while i <= M && j <= M  % Added safeguard for j <= M
        if positions(i) < cumsum_weights(j)
            % Assign index j to resampled particle
            indices(i) = j;
            i = i + 1;  % Move to the next position
        else
            j = j + 1;  % Move to the next cumulative sum
        end
    end
    % Ensure all positions are assigned, handle boundary case
    indices(i:end) = M;  % If j reaches M, assign remaining positions to the last particle
end
% Plot results for SOC estimation using dynamic color generation
num_colors = length(voltage_noise_levels) * length(current_noise_levels);
colors = lines(num_colors);  % Generate distinct colors dynamically
figure;
subplot(3, 1, 1);
plot_idx = 1;
legend_labels = {};  % Initialize empty cell array for legend labels
% Plot Gaussian noise results
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t, squeeze(SOC_est_results_gaussian(v, c, :)), 'Color', colors(plot_idx, :));
        hold on;
        legend_labels{plot_idx} = sprintf('Voltage Noise %.2f, Current Noise %.2f', voltage_noise_levels(v), current_noise_levels(c));  % Create label
        plot_idx = plot_idx + 1;
    end
end
plot(t, SOC_true, 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'True SOC';
xlabel('Time (s)');
ylabel('SOC (%)');
legend(legend_labels, 'Location', 'southwest');
title('SOC Estimation with Gaussian Noise');
grid on;
subplot(3, 1, 2);
plot_idx = 1;
legend_labels = {};  % Reset legend labels
% Plot Pink noise results
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t, squeeze(SOC_est_results_pink(v, c, :)), 'Color', colors(plot_idx, :));
        hold on;
        legend_labels{plot_idx} = sprintf('Voltage Noise %.2f, Current Noise %.2f', voltage_noise_levels(v), current_noise_levels(c));  % Create label
        plot_idx = plot_idx + 1;
    end
end
plot(t, SOC_true, 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'True SOC';
xlabel('Time (s)');
ylabel('SOC (%)');
legend(legend_labels, 'Location', 'southwest');
title('SOC Estimation with Pink Noise');
grid on;
subplot(3, 1, 3);
plot_idx = 1;
legend_labels = {};  % Initialize empty cell array for legend labels
% Plot white noise results
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t, squeeze(SOC_est_results_white(v, c, :)), 'Color', colors(plot_idx, :));
        hold on;
        legend_labels{plot_idx} = sprintf('Voltage Noise %.2f, Current Noise %.2f', voltage_noise_levels(v), current_noise_levels(c));  % Create label
        plot_idx = plot_idx + 1;
    end
end
plot(t, SOC_true, 'k--', 'LineWidth', 2);
legend_labels{end+1} = 'True SOC';
xlabel('Time (s)');
ylabel('SOC (%)');
legend(legend_labels, 'Location', 'southwest');
title('SOC Estimation with White Noise');
grid on;
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
