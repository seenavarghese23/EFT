% Parameters
R0 = 0.01;      % Ohmic resistance
R1 = 0.01;      % First RC pair resistance
C1 = 2000;      % First RC pair capacitance
R2 = 0.01;      % Second RC pair resistance
C2 = 2000;      % Second RC pair capacitance
I = 1; % Current in Amperes

% Current for normal discharge
I_normal = 1; % Current in Amperes

% Increased current for fast discharge
I_high = 10; % Increased current in Amperes

% Time vectors
tspan_normal = [0 10]; % 0 to 10 seconds for normal observation
tspan_fast = [0 1]; % 0 to 1 second for faster observation

% Initial conditions [U1(0), U2(0)]
initial_conditions = [0, 0];

% Differential equations for normal discharge
odefun_normal = @(t, U) [
    (I_normal - U(1)/R1) / C1;
    (I_normal - U(2)/R2) / C2;
];

% Differential equations for high current discharge
odefun_high = @(t, U) [
    (I_high - U(1)/R1) / C1;
    (I_high - U(2)/R2) / C2;
];

% Solve ODE for normal discharge
[t_normal, U_normal] = ode45(odefun_normal, tspan_normal, initial_conditions);
V_true_normal = I_normal * R0 + U_normal(:,1) + U_normal(:,2);

% Solve ODE for fast discharge (shorter time span)
[t_fast, U_fast] = ode45(odefun_normal, tspan_fast, initial_conditions);
V_true_fast_time = I_normal * R0 + U_fast(:,1) + U_fast(:,2);

% Solve ODE for high current discharge
[t_high, U_high] = ode45(odefun_high, tspan_normal, initial_conditions);
V_true_high = I_high * R0 + U_high(:,1) + U_high(:,2);

% Invert the voltage values for plotting
U1_normal_inverted = -U_normal(:,1);
U2_normal_inverted = -U_normal(:,2);
V_true_normal_inverted = -V_true_normal;

U1_fast_inverted = -U_fast(:,1);
U2_fast_inverted = -U_fast(:,2);
V_true_fast_time_inverted = -V_true_fast_time;

U1_high_inverted = -U_high(:,1);
U2_high_inverted = -U_high(:,2);
V_true_high_inverted = -V_true_high;

% Plot comparison results
figure;
subplot(3,1,1);
plot(t_normal, U1_normal_inverted, 'b', t_fast, U1_fast_inverted, 'r', t_high, U1_high_inverted, 'g');
xlabel('Time (s)');
ylabel('Voltage U1 (V)');
title('Voltage across R1');
legend('Normal Discharge', 'Fast Discharge (Short Time)', 'Fast Discharge (High Current)');

subplot(3,1,2);
plot(t_normal, U2_normal_inverted, 'b', t_fast, U2_fast_inverted, 'r', t_high, U2_high_inverted, 'g');
xlabel('Time (s)');
ylabel('Voltage U2 (V)');
title('Voltage across R2');
legend('Normal Discharge', 'Fast Discharge (Short Time)', 'Fast Discharge (High Current)');

subplot(3,1,3);
plot(t_normal, V_true_normal_inverted, 'b', t_fast, V_true_fast_time_inverted, 'r', t_high, V_true_high_inverted, 'g');
xlabel('Time (s)');
ylabel('Terminal Voltage V_{true} (V)');
title('Terminal Voltage V_{true}');
legend('Normal Discharge', 'Fast Discharge (Short Time)', 'Fast Discharge (High Current)');