% Particle Filter for SOC Estimation with 2RC Battery Model including White, Gaussian, and Pink Noise

% Number of particles
num_particles = 1000;

% Initialize particles and weights
particles = repmat([1; 0; 0], 1, num_particles) + 0.1 * randn(3, num_particles);  % Perturb the initial state for diversity
weights = ones(1, num_particles) / num_particles;

% Initialize arrays to store results for different noise types
SOC_est_results_pf_white = zeros(length(voltage_noise_levels), length(current_noise_levels), N);
SOC_est_results_pf_gaussian = zeros(length(voltage_noise_levels), length(current_noise_levels), N);
SOC_est_results_pf_pink = zeros(length(voltage_noise_levels), length(current_noise_levels), N);

% Particle Filter process
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
        V_meas_pink = V_true + voltage_noise_std * generate_pink_noise(N);
        I_meas_pink = I_true + current_noise_std * generate_pink_noise(N);

        % Particle Filter loop for White Noise
        for k = 1:N
            % Propagate particles
            for p = 1:num_particles
                % Predict each particle state using the state transition model
                particles(:, p) = state_transition_jacobian(particles(:, p), I_meas_white(k), dt, Cn, R1, C1, R2, C2);
                
                % Add process noise
                particles(:, p) = particles(:, p) + mvnrnd([0; 0; 0], Q)';
            end
            
            % Update particle weights based on measurement likelihood
            for p = 1:num_particles
                % Predict measurement for each particle
                [z_pred, ~] = measurement_model_jacobian(particles(:, p), I_meas_white(k), R0, R1, R2);
                z_true = [V_meas_white(k); I_meas_white(k)];
                
                % Calculate likelihood
                likelihood = mvnpdf(z_true, z_pred, R_init);
                
                % Update weight
                weights(p) = weights(p) * likelihood;
            end
            
            % Normalize weights
            weights = weights / sum(weights);
            
            % Resample particles if necessary (e.g., using systematic resampling)
            effective_num_particles = 1 / sum(weights.^2);
            if effective_num_particles < num_particles / 2
                resample_idx = resample_particles(weights);
                particles = particles(:, resample_idx);
                weights = ones(1, num_particles) / num_particles;
            end
            
            % Estimate the state (e.g., weighted mean of particles)
            SOC_est_results_pf_white(v, c, k) = sum(particles(1, :) .* weights);
        end

        % Repeat the Particle Filter loop for Gaussian Noise
        for k = 1:N
            % Propagate particles
            for p = 1:num_particles
                % Predict each particle state using the state transition model
                particles(:, p) = state_transition_jacobian(particles(:, p), I_meas_gaussian(k), dt, Cn, R1, C1, R2, C2);
                
                % Add process noise
                particles(:, p) = particles(:, p) + mvnrnd([0; 0; 0], Q)';
            end
            
            % Update particle weights based on measurement likelihood
            for p = 1:num_particles
                % Predict measurement for each particle
                [z_pred, ~] = measurement_model_jacobian(particles(:, p), I_meas_gaussian(k), R0, R1, R2);
                z_true = [V_meas_gaussian(k); I_meas_gaussian(k)];
                
                % Calculate likelihood
                likelihood = mvnpdf(z_true, z_pred, R_init);
                
                % Update weight
                weights(p) = weights(p) * likelihood;
            end
            
            % Normalize weights
            weights = weights / sum(weights);
            
            % Resample particles if necessary (e.g., using systematic resampling)
            effective_num_particles = 1 / sum(weights.^2);
            if effective_num_particles < num_particles / 2
                resample_idx = resample_particles(weights);
                particles = particles(:, resample_idx);
                weights = ones(1, num_particles) / num_particles;
            end
            
            % Estimate the state (e.g., weighted mean of particles)
            SOC_est_results_pf_gaussian(v, c, k) = sum(particles(1, :) .* weights);
        end

        % Repeat the Particle Filter loop for Pink Noise
        for k = 1:N
            % Propagate particles
            for p = 1:num_particles
                % Predict each particle state using the state transition model
                particles(:, p) = state_transition_jacobian(particles(:, p), I_meas_pink(k), dt, Cn, R1, C1, R2, C2);
                
                % Add process noise
                particles(:, p) = particles(:, p) + mvnrnd([0; 0; 0], Q)';
            end
            
            % Update particle weights based on measurement likelihood
            for p = 1:num_particles
                % Predict measurement for each particle
                [z_pred, ~] = measurement_model_jacobian(particles(:, p), I_meas_pink(k), R0, R1, R2);
                z_true = [V_meas_pink(k); I_meas_pink(k)];
                
                % Calculate likelihood
                likelihood = mvnpdf(z_true, z_pred, R_init);
                
                % Update weight
                weights(p) = weights(p) * likelihood;
            end
            
            % Normalize weights
            weights = weights / sum(weights);
            
            % Resample particles if necessary (e.g., using systematic resampling)
            effective_num_particles = 1 / sum(weights.^2);
            if effective_num_particles < num_particles / 2
                resample_idx = resample_particles(weights);
                particles = particles(:, resample_idx);
                weights = ones(1, num_particles) / num_particles;
            end
            
            % Estimate the state (e.g., weighted mean of particles)
            SOC_est_results_pf_pink(v, c, k) = sum(particles(1, :) .* weights);
        end
    end
end

% Plotting SOC Estimation Results for Particle Filter

figure;

% Plot SOC estimation for white noise
subplot(3, 1, 1);
plot(t, SOC_true, 'k', 'LineWidth', 2); hold on;
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        % Plot SOC estimates with white noise
        plot(t, squeeze(SOC_est_results_pf_white(v, c, :)), ...
            'DisplayName', ...
            ['V_{noise} = ', num2str(voltage_noise_levels(v)), ...
             ', I_{noise} = ', num2str(current_noise_levels(c))]);
    end
end
xlabel('Time (s)');
ylabel('SOC');
title('Particle Filter SOC Estimation with White Noise');
legend('True SOC', 'Location', 'Best');
grid on;

% Plot SOC estimation for Gaussian noise
subplot(3, 1, 2);
plot(t, SOC_true, 'k', 'LineWidth', 2); hold on;
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        % Plot SOC estimates with Gaussian noise
        plot(t, squeeze(SOC_est_results_pf_gaussian(v, c, :)), ...
            'DisplayName', ...
            ['V_{noise} = ', num2str(voltage_noise_levels(v)), ...
             ', I_{noise} = ', num2str(current_noise_levels(c))]);
    end
end
xlabel('Time (s)');
ylabel('SOC');
title('Particle Filter SOC Estimation with Gaussian Noise');
legend('True SOC', 'Location', 'Best');
grid on;

% Plot SOC estimation for pink noise
subplot(3, 1, 3);
plot(t, SOC_true, 'k', 'LineWidth', 2); hold on;
for v = 1:length(voltage_noise_levels)
    for c = 1:length(current_noise_levels)
        % Plot SOC estimates with pink noise
        plot(t, squeeze(SOC_est_results_pf_pink(v, c, :)), ...
            'DisplayName', ...
            ['V_{noise} = ', num2str(voltage_noise_levels(v)), ...
             ', I_{noise} = ', num2str(current_noise_levels(c))]);
    end
end
xlabel('Time (s)');
ylabel('SOC');
title('Particle Filter SOC Estimation with Pink Noise');
legend('True SOC', 'Location', 'Best');
grid on;

% Function to generate pink noise
function pink_noise = generate_pink_noise(N)
    % Generate pink noise using Voss-McCartney algorithm
    num_sources = 16;
    rand_sources = rand(num_sources, N);
    pink_noise = sum(rand_sources, 1) - num_sources / 2;
    pink_noise = pink_noise / max(abs(pink_noise));  % Normalize
end

% Function to perform systematic resampling
function resample_idx = resample_particles(weights)
    % Systematic resampling algorithm
    N = length(weights);
    resample_idx = zeros(1, N);
    positions = (rand + (0:(N-1))) / N;
    cumulative_sum = cumsum(weights);
    i = 1;
    j = 1;
    while i <= N
        if positions(i) < cumulative_sum(j)
            resample_idx(i) = j;
            i = i + 1;
        else
            j = j + 1;
        end
    end
end

% State Transition Function for a 2RC Battery Model
function next_state = state_transition_jacobian(current_state, I_meas, dt, Cn, R1, C1, R2, C2)
    % Extract current state variables
    SOC = current_state(1);  % State of charge
    V1 = current_state(2);   % Voltage across the first RC circuit
    V2 = current_state(3);   % Voltage across the second RC circuit

    % State transition equations
    dSOC = -(I_meas * dt) / Cn;  % SOC update (integrating current)
    dV1 = -(V1 / (R1 * C1)) * dt + (I_meas / C1) * dt;  % First RC circuit voltage update
    dV2 = -(V2 / (R2 * C2)) * dt + (I_meas / C2) * dt;  % Second RC circuit voltage update

    % Calculate next state
    SOC_next = SOC + dSOC;
    V1_next = V1 + dV1;
    V2_next = V2 + dV2;

    % Combine into the next state vector
    next_state = [SOC_next; V1_next; V2_next];
end

% Measurement Model Function for the 2RC Battery Model
function [z_pred, H] = measurement_model_jacobian(state, I_meas, R0, R1, R2)
    % Extract state variables
    SOC = state(1);  % State of charge
    V1 = state(2);   % Voltage across the first RC circuit
    V2 = state(3);   % Voltage across the second RC circuit

    % Measurement prediction (voltage)
    V_pred = SOC - I_meas * R0 - V1 - V2;
    z_pred = [V_pred; I_meas];  % Predicted measurement vector

    % Jacobian matrix (partial derivatives of measurement function with respect to state)
    H = [1, -1, -1;
         0,  0,  0];  % Simplified example
end
