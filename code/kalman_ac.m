%% MATLAB Pre-Calculation for FPGA (Corrected with B Matrix)
clear; clc;

FixedPointPrecision = fixdt(1, 32, 20);

% 1. Define Model (From your code)
dt = 0.02;
T = 10;

ampli = 1;

F = single([1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]);
B = single([0.5*dt^2 0; 0 0.5*dt^2; dt 0; 0 dt]); % 4x2 Matrix
H = single(diag([1,1,0,0]));
Q = single(0.1 * eye(4));
R = single(diag([1, 1, ampli, ampli]));

time = single(0:dt:T);

x = [0;0;0;0];
u = single(ones(2, size(time, 2))');

for i = 1:size(time, 2)
    x = [x, F*x(:,i) + B*u(i, :)'];
end
x = x(:, 2:end);
z = x + ampli * randn(4, size(time,2));
z = H*z;
z=z';
