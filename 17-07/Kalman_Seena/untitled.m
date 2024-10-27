% Define state-space matrices
pInitialization.A = [1.1269, -0.4940, 0.1129;
     1.0000, 0, 0;
     0, 1.0000, 0];
pInitialization.B = [-0.3832;
     0.5919;
     0.5191];
pInitialization.C = 1;
pInitialization.D = 0;
pInitialization.Q = 0;
pInitialization.R =[1000, 0, 0;
               0, 10, 0; % Low noise for good V_C1 sensor
               0, 0, 1]; % Low noise for good current sensor;
pInitialization.X0 = [1; 0; 0]; % Initial state (SOC, V_C1, V_C2)
pInitialization.P0 = eye(3); % Initial error covariance
pInitialization.H = 0;
pInitialization.N = 0;
pInitialization.G = 0;
 pInitialization.Z = 0;
 pInitialization.Ts =-1;
 pInitialization.Ns = 3;
TimeDomain = "Continous";
 pInitialization.Nu = 1;
 pInitialization.Nw = 3;
 pInitialization.Ny = 1;
 pInitialization.isSqrtUsed = 0;
 pInitialization.L= [0.78272;0.79185;0.16292];
 pInitialization.M= [0;0;0];

% Control input and measurement
u = [0,0,0]; % Example control input
y= [1,0,0];

% Save variables to workspace
