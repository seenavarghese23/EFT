% Define state-space matrices
A = [1.1269, -0.4940, 0.1129;
     1.0000, 0, 0;
     0, 1.0000, 0];
B = [-0.3832,0,0;
     0.5919,0,0;
     0.5191,0,0];
C = [1;0;0];
D = [0;0;0];
Q = 0;
R = [1,0, 0.1;0,0,0;0,0,0]; % Low noise for good current sensor

X0 = [1; 0; 0]; % Initial state (SOC, V_C1, V_C2)
P0 = eye(3); % Initial error covariance
H = 0;
N = 0;
G = 0;
 Z = 0;
 Ts =-1;
 Ns = 3;
TimeDomain = "Continous";
 Nu = 1;
 Nw = 3;
 Ny = 1;
 isSqrtUsed = 0;
 %L= [0.78272;0.79185;0.16292];
 %M= [0;0;0];

% Control input and measurement
u = [0,0,0;0,0,0;0,0,0]; % Example control input

y= [1,0,0;0,0,0;0,0,0];
% Save variables to workspace
