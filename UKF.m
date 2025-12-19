clear
close all
clc

% Simulation parameters (kept)
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

% Measurement matrix
H = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];

% Covariance matrices
Q = 0.01 * eye(4);
R = diag([measurement_noise_std_gps^2, measurement_noise_std_gps^2, ...
          measurement_noise_std_speed_mps^2, measurement_noise_std_speed_mps^2]);

% True initial state
true_position = [-28; 2];
true_velocity = [speed_mps; speed_mps];
true_positions = [];
measurements = [];
filtered_estimates = [];

% UKF parameters
n = 4;
alpha = 1e-3;
kappa = 0;
beta = 2;
lambda = alpha^2*(n + kappa) - n;
gamma = sqrt(n + lambda);

% Weights
Wm = [lambda/(n+lambda), repmat(1/(2*(n+lambda)),1,2*n)];
Wc = Wm;
Wc(1) = Wc(1) + (1 - alpha^2 + beta);

for k = 1:num_steps
    fprintf('The values for the step %d\n',k);
    disp('---------------------------------');
    % Simulate random changes in velocity
    random_velocity_change = (rand(2, 1) - 0.5) * 2;
    random_acceleration = 0.5 * (randn(2, 1));
    true_velocity = true_velocity + random_velocity_change + random_acceleration * dt;
    disp('Updated velocity is:');
    disp(true_velocity);
    
    % True state (nonlinear)
    true_position(1) = true_position(1) + true_velocity(1)*dt + 0.1*sin(true_position(1));
    true_position(2) = true_position(2) + true_velocity(2)*dt + 0.1*sin(true_position(2));
    disp('Updated position is:');
    disp(true_position);
    
    % Measurements
    gps_measurement = true_position + measurement_noise_std_gps * randn(2, 1);
    disp('Random GPS noise:');
    disp(gps_measurement);
    speedometer_measurement = true_velocity + measurement_noise_std_speed_mps * randn(2, 1);
    disp('Random speed noise');
    disp(speedometer_measurement);
    
    measurements = [measurements; gps_measurement', speedometer_measurement'];
    
    % ----------------- UKF PREDICTION -----------------
    % sigma points generation
    S = chol(P, 'lower');
    sigma_pts = zeros(n, n+2);
    sigma_pts(:,1) = x_est;
    for i = 1:n
        sigma_pts(:,1+i)   = x_est + gamma * S(:,i);
        sigma_pts(:,1+n+i) = x_est - gamma * S(:,i);
    end
    
    % propagate sigma points through nonlinear motion model f
    sigma_pred = zeros(n, 2*n+1);
    for i = 1:(2*n+1)
        x = sigma_pts(:,i);
        % f(x): same nonlinearity as used for true state
        xpred = [ x(1) + x(3)*dt + 0.1*sin(x(1));
                  x(2) + x(4)*dt + 0.1*sin(x(2));
                  x(3);
                  x(4) ];
        sigma_pred(:,i) = xpred; % no control (u==0)
    end
    
    % predicted mean
    x_pred = zeros(n,1);
    for i = 1:(2*n+1)
        x_pred = x_pred + Wm(i) * sigma_pred(:,i);
    end
    
    % predicted covariance
    P_pred = Q;
    for i = 1:(2*n+1)
        d = sigma_pred(:,i) - x_pred;
        P_pred = P_pred + Wc(i) * (d * d');
    end
    disp('The next state (predicted) is:');
    disp(x_pred);
    disp('The initial P matrix (predicted) is:');
    disp(P_pred);
    
    % ----------------- UKF UPDATE -----------------
    % measurement sigma points (measurement function is identity here)
    z_sigma = sigma_pred; % because measurement is direct state (H = I)
    
    % predicted measurement
    z_pred = zeros(n,1);
    for i = 1:(2*n+1)
        z_pred = z_pred + Wm(i) * z_sigma(:,i);
    end
    
    % innovation covariance and cross-cov
    Pzz = R;
    Pxz = zeros(n);
    for i = 1:(2*n+1)
        dz = z_sigma(:,i) - z_pred;
        dx = sigma_pred(:,i) - x_pred;
        Pzz = Pzz + Wc(i) * (dz * dz');
        Pxz = Pxz + Wc(i) * (dx * dz');
    end
    
    % Kalman gain
    K = Pxz / Pzz;
    disp('Kalman gain is:');
    disp(K);
    
    measurement = [gps_measurement; speedometer_measurement];
    x_est = x_pred + K * (measurement - z_pred);
    disp('The estimated state is:');
    disp(x_est);
    P = P_pred - K * Pzz * K';
    disp('The P covariance matrix is:');
    disp(P);
    disp('---------------------------------');
    
    % Store
    true_positions = [true_positions; true_position'];
    filtered_estimates = [filtered_estimates; x_est(1:4)'];
end

% Plot results (matching style)
time = (1:num_steps) * dt;

figure;
subplot(2, 1, 1);
plot(true_positions(:, 1), true_positions(:, 2), '-g', 'LineWidth', 1.5); hold on;
plot(measurements(:, 1), measurements(:, 2), 'xr'); hold on;
plot(filtered_estimates(:, 1), filtered_estimates(:, 2), '-b', 'LineWidth', 1.5);
legend('True Position', 'GPS Measurements', 'UKF Estimated Position');
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('UKF: 2D Position Estimation (Random Motion, small nonlinearity)');
grid on;

subplot(2, 1, 2);
vx_filtered = filtered_estimates(:,3);
plot(time, vx_filtered, '-b', 'LineWidth', 1.5); hold on;
plot(time, measurements(:, 3), 'xr');
legend('UKF Filtered Vx', 'Speedometer Measurements (Vx)');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('UKF: Velocity Estimation (Vx)');
grid on;
