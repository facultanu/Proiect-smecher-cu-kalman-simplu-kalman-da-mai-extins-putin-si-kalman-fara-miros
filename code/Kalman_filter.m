% Kalman Filter Example: Car moving randomly in 2D

clear
close all
clc

% Simulation parameters
dt = 0.2;                % Time step (seconds)
num_steps = 60/dt;        % 1 minute of motion (60 seconds)
speed_kmh = 80;        % Initial Speed in km/h (can be modified)
speed_mps = speed_kmh * 1000 / 3600; % Initial Speed in m/s (16.67 m/s)
measurement_noise_std_gps = 1;  % GPS measurement noise (meters)
measurement_noise_std_speed = 2; % Speedometer noise (km/h)
measurement_noise_std_speed_mps = measurement_noise_std_speed * 1000 / 3600; % Speedometer noise (m/s)

% Initial state estimate
x_est = [-20; -10; 0; 0];  % Initial state: [x; y; vx; vy]
P = eye(4);            % Initial covariance matrix (uncertainty about initial state)

% State transition matrix (models motion with constant velocity)
A = [1 0 dt 0;  % Position x depends on velocity x
     0 1 0 dt;  % Position y depends on velocity y
     0 0 1 0;  % Velocity x remains constant
     0 0 0 1]; % Velocity y remains constant

% Control input model
B = [0.5*dt^2 0;
    0 0.5*dt^2;
        dt 0;
        0 dt]; % Acceleration (x, y) control matrix
u = [0; 0]; % Initial acceleration (no control at start)

% Measurement matrix (we measure position and velocity separately)
H = [1 0 0 0;   % GPS measures x position
     0 1 0 0;   % GPS measures y position
     0 0 1 0;   % Speedometer measures x velocity
     0 0 0 1];  % Speedometer measures y velocity

% Covariance matrices
Q = 0.01 * eye(4); % Process noise covariance (uncertainty in motion model)
R = diag([measurement_noise_std_gps^2, measurement_noise_std_gps^2, ...
          measurement_noise_std_speed_mps^2, measurement_noise_std_speed_mps^2]); % Measurement noise covariance

%% Simulate true position and velocity
true_position = [-28; 2];                 % Start position in x and y
true_velocity = [speed_mps; speed_mps]; % Initial velocity (can be randomized)
true_positions = [];
measurements = [];                      % Noisy GPS and speedometer readings
filtered_estimates = [];                % Kalman Filter's estimate of position and velocity

for k = 1:num_steps
    fprintf('The values for the step %d\n',k);
    disp('---------------------------------');
    % Simulate random changes in velocity
    random_velocity_change = (rand(2, 1) - 0.5) * 2;  % Random change in velocity
    random_acceleration = 0.5 * (randn(2, 1)); % Random acceleration (m/s^2)
    true_velocity = true_velocity + random_velocity_change + random_acceleration * dt; % Update velocity with random acceleration
    disp('Updated velocity is:');
    disp(true_velocity);
    
    % Simulate true position and velocity
    true_position = true_position + true_velocity * dt; % Update true position
    disp('Updated position is:');
    disp(true_position);
    
    % Simulate noisy measurements
    gps_measurement = true_position + measurement_noise_std_gps * randn(2, 1); % Noisy GPS for x and y
    disp('Random GPS noise:');
    disp(gps_measurement);
    speedometer_measurement = true_velocity + measurement_noise_std_speed_mps * randn(2, 1); % Noisy speedometer for x and y
    disp('Random speed noise');
    disp(speedometer_measurement);
    
    % Store measurements
    measurements = [measurements; gps_measurement', speedometer_measurement'];
    
    % Prediction step
    x_pred = A * x_est + B * u;
    disp('The next state (predicted) is:');
    disp(x_pred);
    P_pred = A * P * A' + Q;
    disp('The initial P matrix (predicted) is:');
    disp(P_pred);
    
    % Update step
    measurement = [gps_measurement; speedometer_measurement]; % Measurement vector
    K = P_pred * H' / (H * P_pred * H' + R); % Compute Kalman gain
    disp('Kalman gain is:');
    disp(K);
    x_est = x_pred + K * (measurement - H * x_pred); % Update state estimate
    disp('The estimated state is:');
    disp(x_est);
    P = (eye(4) - K * H) * P_pred;      % Update covariance matrix
    disp('The P covariance matrix is:');
    disp(P);

    disp('---------------------------------');
    
    % Store true and filtered positions
    true_positions = [true_positions; true_position'];
    filtered_estimates = [filtered_estimates; x_est(1:4)'];
end

% Plot results
time = (1:num_steps) * dt; % Time in seconds

% Position plot (x and y)
figure;
subplot(2, 1, 1);
plot(true_positions(:, 1), true_positions(:, 2), '-g', 'LineWidth', 1.5); hold on;
plot(measurements(:, 1), measurements(:, 2), 'xr'); hold on;
plot(filtered_estimates(:, 1), filtered_estimates(:, 2), '-b', 'LineWidth', 1.5);
legend('True Position', 'GPS Measurements', 'Filtered Position');
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Kalman Filter: 2D Position Estimation (Random Motion)');
grid on;

% Velocity plot (x and y)
subplot(2, 1, 2);
plot(time, true_velocity(1) * ones(num_steps, 1), '-g', 'LineWidth', 1.5); hold on;
plot(time, measurements(:, 3), 'xr'); hold on;
plot(time, filtered_estimates(:, 3), '-b', 'LineWidth', 1.5);
legend('True Velocity', 'Speedometer Measurements', 'Filtered Velocity');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Kalman Filter: 2D Velocity Estimation (Random Motion)');
grid on;
