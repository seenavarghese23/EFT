% Define the measurement noise standard deviations
sigma_V = 0.01; % Standard deviation for voltage sensor
sigma_I = 0.005; % Standard deviation for current sensor

% Construct the measurement noise covariance matrix R
R = [sigma_V^2, 0;
     0, sigma_I^2];
Q=1;
Q_2RC=1;
Q2RC=eye(size(Q_2RC,1));

%% FC parameters
Fuel_Mass = 2*15.51*1000; % Fuel Mass in kg

%% Cell parameters
%Note: set polynomial coefficents (in the OCV mathfunction block of the model) correlating of the loaded cell
%load Cell_data\LG_INR_18650_MJ1.mat;
%load Cell_data\Sony_Konion_US18650_VTC5A.mat;
%load Cell_data\Murata_US1_18650_VTC6.mat;
addpath 'C:\Users\DELL\Downloads\SoC_Est'\Cell_data\
load Molicell_INR_21700_P45B.mat;

C1 = 0.9082;
C2 = 1.1322;
int_cell_resistance = 0.0090;
R1 = 0.001;
R2 = 0.001;

stepSize=0.5;
%% Fuel Cell Parameters
P_op = 1.5; % Operating Temp. of Fuel Cell in bar
T_op = 65; % Operating Temp of Fuel Cell in Celcius
n_cell_stack = 500; % number of cells per stack
% *Include in model and a variable -- Tank with fuel consumption and initial fuel variable input

%% General cell/battery parameters
Module_count = 8; % Left for user to adjust accordingly for specific number of modules
series_cells = Module_count * 16;
parallel_cells = 11; %Number of modules multiplied by number of parallel cells
init_cell_temp = 25 + 273;
nominal_ambient_temp = 25 + 273;
electronic_resistance = 0.00001;
initial_state_of_charge = 90; %initial state of charge as percent
weight_integration_factor = 1.5;
cell_cp_heat = 950; % J/ kg K
cell_area = 2 * pi * 9e-3 * 65e-3; % = 2 * pi * r * h
h_conv = 40; % W / m^2 K % 6 bis 8 for natural convection, greater than 10 for forced cooling

%% General pack calculations
operating_voltage_pack = nom_cell_voltage * series_cells;
maximum_capacity_pack =  max_cell_capacity * parallel_cells;
cut_off_voltage_pack = cutoff_cell_voltage * series_cells;
fully_charged_voltage_pack = fullycha_cell_voltage * series_cells;
rated_capacity_pack =  rated_cell_capacity * parallel_cells;
rated_capacity_pack_105 = rated_capacity_pack * 1.05;
internal_resistance_pack = int_cell_resistance * series_cells / parallel_cells + electronic_resistance;
weight_pack = weight_cell * series_cells * parallel_cells * weight_integration_factor;
cell_area_pack = cell_area * series_cells * parallel_cells;

%% 2RC-Thevenin Model calculations
R1_pack = R1 * series_cells / parallel_cells;   
C1_pack = C1 * parallel_cells / series_cells;
R2_pack = R2 * series_cells / parallel_cells;
C2_pack = C2 * parallel_cells / series_cells;




