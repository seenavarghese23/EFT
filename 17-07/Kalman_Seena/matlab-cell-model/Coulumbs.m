% Simulate the model and retrieve the output
simOut = sim('Coulombs.slx');

% Retrieve the SoC signals from the simulation output
% Assuming the signals are logged as 'SoC' and 'ErroredSoC'
SoC = simOut.get('DatawithoutNoise'); % For 'To Workspace' block with 'SoC' variable name
ErroredSoC = simOut.get('DataNoise'); % For 'To Workspace' block with 'ErroredSoC' variable name

% Display the element names inthe SoC structure
SoC_elements = SoC.getElementNames();
disp(SoC_elements);

% Access the specific elements by name or index
% Example: accessing the first and second signals
signal1 = SoC.getElement(1).Values; % SoC
signal2 = ErroredSoC.getElement(1).Values; % Errored SoC

% Extract time and data
time = signal1.Time; % Assuming both signals have the same time vector
signal1_data = signal1.Data;
signal2_data = signal2.Data;

% Plot the signals to verify
figure;

% Plot the first signal (SoC)
subplot(2,1,1);
plot(time, signal1_data);
title('Signal 1: SoC');
xlabel('Time (s)');
ylabel('SoC');

% Plot the second signal (Errored SoC) along with the first signal for comparison
subplot(2,1,2);
plot(time, signal2_data);
hold on; % Allow multiple plots on the same axes
plot(time, signal1_data);
title('Signal 2: Errored SoC and SoC');
xlabel('Time (s)');
ylabel('SoC');
legend('Errored SoC', 'SoC'); % Add a legend to distinguish the signals
hold off; % Release the hold on the axes
