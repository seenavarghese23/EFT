% UKF for SOC Estimation with 2RC Battery Model including White, Gaussian, and Pink Noise

% Battery parameters
Cn = 4.5 * 3600;  % Nominal capacity in As (Ampere-seconds)
R0 = 0.002;        % Ohmic resistance
R1 = 0.003;        % First RC pair resistance
C1 = 1000e-3;        % First RC pair capacitance
R2 = 0.00522;        % Second RC pair resistance
C2 = 1000e-3;        % Second RC pair capacitance

% Simulation parameters
t_end = 3600;     % Simulation time in seconds
dt = 1;           % Time step in seconds
t = 0:dt:t_end;   % Time vector
N = length(t);    % Number of simulation steps

% Generate simulated current profile (e.g., constant discharge)
I_true = -5 * ones(1, N);  % 1A discharge current

% Generate true SOC profile
%SOC_true = linspace(1, 0.7, N);  % Linear discharge from 100% to 70%
% Load SOC time data
addpath 'C:\Users\DELL\Documents\Master Thesis EFT\github\Old as 17-07-2024\Test'
load SOC_time_data1234.mat;
SOC_true = data1;

% Generate true voltage measurements based on the battery model
V_true = ocv(SOC_true) - I_true*R0 - I_true*R1.*(1-exp(-t./(R1*C1))) - I_true*R2.*(1-exp(-t./(R2*C2)));

% Define noise levels for comparison
voltage_noise_levels = [0.01, 0.05, 0.1];  % Standard deviations for voltage noise
current_noise_levels = [0.001, 0.05, 0.1];  % Standard deviations for current noise

% Initialize arrays to store results for different noise types
SOC_est_results_white = zeros(length(voltage_noise_levels), length(current_noise_levels), N);
SOC_est_results_gaussian = zeros(length(voltage_noise_levels), length(current_noise_levels), N);
SOC_est_results_pink = zeros(length(voltage_noise_levels), length(current_noise_levels), N);

% Helper function to generate pink noise
function x = pinknoise(N)
    white_noise = randn(1, N);
    x = filter([0.021 0.044 0.021], [1 -2.347 1.977 -0.626], white_noise);
end

% UKF Parameters
alpha = 1e-3;       % Default, tunable parameter (typically small)
kappa = 0;          % Secondary scaling parameter (usually 0 or 3-n)
beta = 2;           % Optimal for Gaussian distributions
L = 3;              % Dimension of the state vector

lambda = alpha^2 * (L + kappa) - L;
gamma = sqrt(L + lambda);

% Weights for mean and covariance
Wm = [lambda / (L + lambda), 0.5 / (L + lambda) * ones(1, 2*L)];
Wc = Wm;
Wc(1) = Wc(1) + (1 - alpha^2 + beta);

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

        % Initialize UKF for White Noise
        x_ukf = [100; 0; 0];  % Initial state [SOC; V_RC1; V_RC2]
        P_ukf = diag([0.01^1, 0.01^1, 0.01^1]);  % Initial state covariance
        Q = diag([10e-7, 10e-5, 10e-5]);         % Process noise covariance
        R_ukf = diag([max(voltage_noise_std^2, 1e-10), max(current_noise_std^2, 1e-10)]);  % Measurement noise covariance

        % UKF loop for White Noise
        for k = 1:N
            % Generate sigma points
            X = generate_sigma_points(x_ukf, P_ukf, gamma);
            
            % Predict sigma points
            X_pred = zeros(size(X));
            for i = 1:(2*L + 1)
                X_pred(:, i) = state_transition(X(:, i), I_meas_white(k), dt, Cn, R1, C1, R2, C2);
            end
            
            % Predicted state mean
            x_pred = X_pred * Wm';
            
            % Predicted state covariance
            P_pred = Q;
            for i = 1:(2*L + 1)
                diff = X_pred(:, i) - x_pred;
                P_pred = P_pred + Wc(i) * (diff * diff');
            end
            
            % Transform sigma points to measurement space
            Z_pred = zeros(2, 2*L + 1); % Measurement has two dimensions (voltage and current)
            for i = 1:(2*L + 1)
                Z_pred(:, i) = measurement_model(X_pred(:, i), I_meas_white(k), R0, R1, R2);
            end
            
            % Predicted measurement mean
            z_pred = Z_pred * Wm';
            
            % Innovation covariance
            S = R_ukf;
            for i = 1:(2*L + 1)
                diff = Z_pred(:, i) - z_pred;
                S = S + Wc(i) * (diff * diff');
            end
            
            % Cross-covariance
            Pxz = zeros(L, 2);
            for i = 1:(2*L + 1)
                Pxz = Pxz + Wc(i) * (X_pred(:, i) - x_pred) * (Z_pred(:, i) - z_pred)';
            end
            
            % Kalman gain
            K = Pxz / S;
            
            % Update state estimate
            z = [V_meas_white(k); I_meas_white(k)];
            x_ukf = x_pred + K * (z - z_pred);
            
            % Update state covariance
            P_ukf = P_pred - K * S * K';
            
            % Store results
            SOC_est_results_white(v, c, k) = x_ukf(1);
        end

        % Repeat the UKF loop for Gaussian Noise
        % Initialize UKF for Gaussian Noise
        x_ukf = [100; 0; 0];  % Initial state [SOC; V_RC1; V_RC2]
        P_ukf = diag([0.01^1, 0.01^1, 0.01^1]);  % Initial state covariance
        R_ukf = diag([max(voltage_noise_std^2, 1e-10), max(current_noise_std^2, 1e-10)]);  % Measurement noise covariance

        for k = 1:N
            % Generate sigma points
            X = generate_sigma_points(x_ukf, P_ukf, gamma);
            
            % Predict sigma points
            X_pred = zeros(size(X));
            for i = 1:(2*L + 1)
                X_pred(:, i) = state_transition(X(:, i), I_meas_gaussian(k), dt, Cn, R1, C1, R2, C2);
            end
            
            % Predicted state mean
            x_pred = X_pred * Wm';
            
            % Predicted state covariance
            P_pred = Q;
            for i = 1:(2*L + 1)
                diff = X_pred(:, i) - x_pred;
                P_pred = P_pred + Wc(i) * (diff * diff');
            end
            
            % Transform sigma points to measurement space
            Z_pred = zeros(2, 2*L + 1); % Measurement has two dimensions (voltage and current)
            for i = 1:(2*L + 1)
                Z_pred(:, i) = measurement_model(X_pred(:, i), I_meas_gaussian(k), R0, R1, R2);
            end
            
            % Predicted measurement mean
            z_pred = Z_pred * Wm';
            
            % Innovation covariance
            S = R_ukf;
            for i = 1:(2*L + 1)
                diff = Z_pred(:, i) - z_pred;
                S = S + Wc(i) * (diff * diff');
            end
            
            % Cross-covariance
            Pxz = zeros(L, 2);
            for i = 1:(2*L + 1)
                Pxz = Pxz + Wc(i) * (X_pred(:, i) - x_pred) * (Z_pred(:, i) - z_pred)';
            end
            
            % Kalman gain
            K = Pxz / S;
            
            % Update state estimate
            z = [V_meas_gaussian(k); I_meas_gaussian(k)];
            x_ukf = x_pred + K * (z - z_pred);
            
            % Update state covariance
            P_ukf = P_pred - K * S * K';
            
            % Store results
            SOC_est_results_gaussian(v, c, k) = x_ukf(1);
        end

        % Repeat the UKF loop for Pink Noise
        % Initialize UKF for Pink Noise
        x_ukf = [100; 0; 0];  % Initial state [SOC; V_RC1; V_RC2]
        P_ukf = diag([0.01^1, 0.01^1, 0.01^1]);  % Initial state covariance
        R_ukf = diag([max(voltage_noise_std^2, 1e-10), max(current_noise_std^2, 1e-10)]);  % Measurement noise covariance

        for k = 1:N
            % Generate sigma points
            X = generate_sigma_points(x_ukf, P_ukf, gamma);
            
            % Predict sigma points
            X_pred = zeros(size(X));
            for i = 1:(2*L + 1)
                X_pred(:, i) = state_transition(X(:, i), I_meas_pink(k), dt, Cn, R1, C1, R2, C2);
            end
            
            % Predicted state mean
            x_pred = X_pred * Wm';
            
            % Predicted state covariance
            P_pred = Q;
            for i = 1:(2*L + 1)
                diff = X_pred(:, i) - x_pred;
                P_pred = P_pred + Wc(i) * (diff * diff');
            end
            
            % Transform sigma points to measurement space
            Z_pred = zeros(2, 2*L + 1); % Measurement has two dimensions (voltage and current)
            for i = 1:(2*L + 1)
                Z_pred(:, i) = measurement_model(X_pred(:, i), I_meas_pink(k), R0, R1, R2);
            end
            
            % Predicted measurement mean
            z_pred = Z_pred * Wm';
            
            % Innovation covariance
            S = R_ukf;
            for i = 1:(2*L + 1)
                diff = Z_pred(:, i) - z_pred;
                S = S + Wc(i) * (diff * diff');
            end
            
            % Cross-covariance
            Pxz = zeros(L, 2);
            for i = 1:(2*L + 1)
                Pxz = Pxz + Wc(i) * (X_pred(:, i) - x_pred) * (Z_pred(:, i) - z_pred)';
            end
            
            % Kalman gain
            K = Pxz / S;
            
            % Update state estimate
            z = [V_meas_pink(k); I_meas_pink(k)];
            x_ukf = x_pred + K * (z - z_pred);
            
            % Update state covariance
            P_ukf = P_pred - K * S * K';
            
            % Store results
            SOC_est_results_pink(v, c, k) = x_ukf(1);
        end
    end
end

% Plot SOC estimation results for different noise types
figure;
plot(t, SOC_true, 'k--', 'LineWidth', 2); hold on;
plot(t, squeeze(SOC_est_results_white(1, 1, :)), 'r', 'LineWidth', 1.5);
plot(t, squeeze(SOC_est_results_gaussian(1, 1, :)), 'b', 'LineWidth', 1.5);
plot(t, squeeze(SOC_est_results_pink(1, 1, :)), 'g', 'LineWidth', 1.5);
legend('True SOC', 'White Noise', 'Gaussian Noise', 'Pink Noise');
xlabel('Time (s)');
ylabel('State of Charge (SOC)');
title('SOC Estimation Using UKF for Different Noise Types');
grid on;

% Helper Functions
function X = generate_sigma_points(x, P, gamma)
    L = numel(x); % Number of states
    X = zeros(L, 2 * L + 1); % Sigma points matrix
    X(:, 1) = x; % First column is the mean
    sqrtP = chol(P, 'lower');
    for i = 1:L
        X(:, i + 1) = x + gamma * sqrtP(:, i);
        X(:, i + 1 + L) = x - gamma * sqrtP(:, i);
    end
end

function x_next = state_transition(x, I, dt, Cn, R1, C1, R2, C2)
    SOC = x(1);
    V_RC1 = x(2);
    V_RC2 = x(3);

    dV_RC1 = (-V_RC1 + I * R1) / (R1 * C1);
    dV_RC2 = (-V_RC2 + I * R2) / (R2 * C2);
    
    SOC_next = SOC - I * dt / Cn;
    V_RC1_next = V_RC1 + dV_RC1 * dt;
    V_RC2_next = V_RC2 + dV_RC2 * dt;

    x_next = [SOC_next; V_RC1_next; V_RC2_next];
end

function z = measurement_model(x, I, R0, R1, R2)
    SOC = x(1);
    V_RC1 = x(2);
    V_RC2 = x(3);

    OCV = ocv(SOC);
    V = OCV - I * R0 - V_RC1 - V_RC2;
    
    z = [V; I];
end

function ocv_val = ocv(SOC)
    % Open-circuit voltage (OCV) function
    ocv_val = 3.6 + 0.4 * SOC;  % Example linear relationship
end
% Plotting SOC Estimation Results

% Create a new figure
figure;

% Plot SOC estimation for white noise
subplot(3, 1, 1);
plot(t, SOC_true, 'k', 'LineWidth', 2); hold on;
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t, squeeze(SOC_est_results_white(v, c, :)), 'DisplayName', ...
            ['V_{noise} = ', num2str(voltage_noise_levels(v)), ...
             ', I_{noise} = ', num2str(current_noise_levels(c))]);
    end
end
xlabel('Time (s)');
ylabel('SOC');
title('SOC Estimation with White Noise');
legend('True SOC', 'Location', 'Best');
grid on;

% Plot SOC estimation for Gaussian noise
subplot(3, 1, 2);
plot(t, SOC_true, 'k', 'LineWidth', 2); hold on;
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t, squeeze(SOC_est_results_gaussian(v, c, :)), 'DisplayName', ...
            ['V_{noise} = ', num2str(voltage_noise_levels(v)), ...
             ', I_{noise} = ', num2str(current_noise_levels(c))]);
    end
end
xlabel('Time (s)');
ylabel('SOC');
title('SOC Estimation with Gaussian Noise');
legend('True SOC', 'Location', 'Best');
grid on;

% Plot SOC estimation for pink noise
subplot(3, 1, 3);
plot(t, SOC_true, 'k', 'LineWidth', 2); hold on;
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        plot(t, squeeze(SOC_est_results_pink(v, c, :)), 'DisplayName', ...
            ['V_{noise} = ', num2str(voltage_noise_levels(v)), ...
             ', I_{noise} = ', num2str(current_noise_levels(c))]);
    end
end
xlabel('Time (s)');
ylabel('SOC');
title('SOC Estimation with Pink Noise');
legend('True SOC', 'Location', 'Best');
grid on;

% Adjust layout
sgtitle('UKF SOC Estimation Results for Different Noise Types');
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
        deviation_gaussian(v, c) = sqrt((mean(abs_diff_gaussian).^2));
        
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
