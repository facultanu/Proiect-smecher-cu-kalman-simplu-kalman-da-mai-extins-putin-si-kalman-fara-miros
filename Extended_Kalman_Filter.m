clear
close all
clc

% Simulation parameters (kept exactly as your original)
dt = 0.2;                
num_steps = 60/dt;        
speed_kmh = 80;        
speed_mps = speed_kmh * 1000 / 3600; 
measurement_noise_std_gps = 1;  
measurement_noise_std_speed = 2; 
measurement_noise_std_speed_mps = measurement_noise_std_speed * 1000 / 3600;

% Initial state estimate
x_est = [-20; -10; 0; 0];  
P = eye(4);            

% State transition (used only as a template; EKF uses nonlinear f and Jacobian)
A = [1 0 dt 0;
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];

% Control input model (kept)
B = [0.5*dt^2 0;
     0 0.5*dt^2;
     dt 0;
     0 dt];
u = [0; 0];

% Measurement matrix (same as yours)
H = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

% Covariance matrices
Q = 0.01 * eye(4);
R = diag([measurement_noise_std_gps^2, measurement_noise_std_gps^2, ...
          measurement_noise_std_speed_mps^2, measurement_noise_std_speed_mps^2]);

% True initial state (close to your original)
true_position = [-28; 2];
true_velocity = [speed_mps; speed_mps];
true_positions = [];
measurements = [];
filtered_estimates = [];

for k = 1:num_steps
    fprintf('The values for the step %d\n',k);
    disp('---------------------------------');
    % Simulate random changes in velocity
    random_velocity_change = (rand(2, 1) - 0.5) * 2;
    random_acceleration = 0.5 * (randn(2, 1));
    true_velocity = true_velocity + random_velocity_change + random_acceleration * dt;
    disp('Updated velocity is:');
    disp(true_velocity);
    
    % ----------------- TRUE STATE (nonlinear) -----------------
    % small nonlinearity added to position update
    true_position(1) = true_position(1) + true_velocity(1)*dt + 0.1*sin(true_position(1));
    true_position(2) = true_position(2) + true_velocity(2)*dt + 0.1*sin(true_position(2));
    disp('Updated position is:');
    disp(true_position);
    
    % Simulate noisy measurements
    gps_measurement = true_position + measurement_noise_std_gps * randn(2, 1);
    disp('Random GPS noise:');
    disp(gps_measurement);
    speedometer_measurement = true_velocity + measurement_noise_std_speed_mps * randn(2, 1);
    disp('Random speed noise');
    disp(speedometer_measurement);
    
    % Store measurements
    measurements = [measurements; gps_measurement', speedometer_measurement'];
    
    % ----------------- EKF PREDICTION -----------------
    % Nonlinear motion model f(x)
    fx = [ x_est(1) + x_est(3)*dt + 0.1*sin(x_est(1));
           x_est(2) + x_est(4)*dt + 0.1*sin(x_est(2));
           x_est(3);
           x_est(4) ];
    x_pred = fx + B*u;   % include control (u is zero here)
    disp('The next state (predicted) is:');
    disp(x_pred);
    
    % Jacobian F = df/dx (4x4)
    % df1/dx1 = 1 + 0.1*cos(x1), df1/dvx = dt, similarly for y
    F = [ 1 + 0.1*cos(x_est(1))  0  dt  0;
          0  1 + 0.1*cos(x_est(2))  0  dt;
          0  0  1  0;
          0  0  0  1 ];
    P_pred = F * P * F' + Q;
    disp('The initial P matrix (predicted) is:');
    disp(P_pred);
    
    % ----------------- EKF UPDATE -----------------
    measurement = [gps_measurement; speedometer_measurement];
    K = P_pred * H' / (H * P_pred * H' + R);
    disp('Kalman gain is:');
    disp(K);
    x_est = x_pred + K * (measurement - H * x_pred);
    disp('The estimated state is:');
    disp(x_est);
    P = (eye(4) - K * H) * P_pred;
    disp('The P covariance matrix is:');
    disp(P);
    disp('---------------------------------');
    
    % Store
    true_positions = [true_positions; true_position'];
    filtered_estimates = [filtered_estimates; x_est(1:4)'];
end

% Plot results
time = (1:num_steps) * dt;

figure;
subplot(2, 1, 1);
plot(true_positions(:, 1), true_positions(:, 2), '-g', 'LineWidth', 1.5); hold on;
plot(measurements(:, 1), measurements(:, 2), 'xr'); hold on;
plot(filtered_estimates(:, 1), filtered_estimates(:, 2), '-b', 'LineWidth', 1.5);
legend('True Position', 'GPS Measurements', 'EKF Estimated Position');
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('EKF: 2D Position Estimation (Random Motion, small nonlinearity)');
grid on;

subplot(2, 1, 2);

vx_filtered = filtered_estimates(:,3);
plot(time, vx_filtered, '-b', 'LineWidth', 1.5); hold on;
plot(time, measurements(:, 3), 'xr');
legend('EKF Filtered Vx', 'Speedometer Measurements (Vx)');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('EKF: Velocity Estimation (Vx)');
grid on;
