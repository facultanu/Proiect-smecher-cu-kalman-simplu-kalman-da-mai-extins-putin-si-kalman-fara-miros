%% MATLAB Pre-Calculation for FPGA (Corrected with B Matrix)
clear; clc;

% 1. Define Model (From your code)
dt = 0.02;
T = 10;

F = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
B = [0.5*dt^2 0; 0 0.5*dt^2; dt 0; 0 dt]; % 4x2 Matrix
H = eye(4);
Q = 0.01 * eye(4);
R = diag([1, 1, 0.55, 0.55]);

time = 0:dt:T;
z = [1;1;1;1]*sin(time) + 0.1*randn(4, size(time, 2));
z=z';
u = zeros(2, size(time, 2))';