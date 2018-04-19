% System Model
%
%   dx(t)/dt = a(t)*x(t) + b(t)*u(t) + d(t),
%   x(0) = 10;
%   da(t)/dt = 0.05*sin(3t) + w1(t),
%   a(0) = 2.5;
%   b(t) = 2;
%   d(t) = 1.5*sin(t) + 1 + 2*w2(t);
%   where w1(t), w2(t) ~ N(0,1);
%

%% Clear and Start Clocking
clc; clear; close all; tic;

%% Preprare Simulation

% the differential step size
dt = 0.0001; N = 120000;

% preallocation
xref = zeros(1,N+2);
u_compensate = zeros(1,N+2); u_pid = zeros(1,N+2); u_all = zeros(1,N+2);
a = zeros(1,N+2); b = zeros(1,N+2); d = zeros(1,N+2);
dx_dt = zeros(1,N+2); x = zeros(1,N+2);
RLS_output = zeros(1,N+2); RLS_regressor = zeros(1,N+2); 
RLS_error = zeros(1,N+2);
RLS_g = zeros(1,N+2); RLS_P = zeros(1,N+2); 
RLS_hata = zeros(1,N+2);
ESO_error = zeros(1,N+2); 
ESO_hatdx_dt = zeros(1,N+2); ESO_hatdd_dt = zeros(1,N+2);
ESO_hatx = zeros(1,N+2); ESO_hatd = zeros(1,N+2);

% the true parameters
a(1) = 2.5; b_N = 2;

% the tunable parameters of RLS estimator
lambda = 0.985; RLS_P(1) = 100; 
RLS_hata(1) = 0;

% the tunable parameters of ESO
pole_tuning = 2; k1 = 2*pole_tuning; k2 = pole_tuning^2;
ESO_hatx(1) = 0; ESO_hatd(1) = 0;

% the tunable parameters of PID
kp = 2; ki = 0; kd = 0;

% the initial value
x(1) = 10;

for k = 1:N+1
    
    %% Reference
    xref(k) = 0;
    
    %% Controller
    u_compensate(k) = ( -RLS_hata(k)*x(k)-ESO_hatd(k) )/b_N;
    u_pid(k) = ( kp*(xref(k)-x(k)) )/b_N;
    u_all(k) = u_pid(k)+u_compensate(k);
    
    %% System
    % time-varying parameter
    a(k+1) = a(k) + dt*( 1*randn+0.5*sin(3*dt*k) );
    % disturbance
    d(k) = 2.5*sin(5*dt*k) + 4 + 1*randn;
    % dynamics
    dx_dt(k) = a(k)*x(k) + b_N*u_all(k) + d(k);
    x(k+1) = x(k) + dt*dx_dt(k);
    
    %% RLS
    RLS_output(k) = dx_dt(k) - b_N*u_all(k) - ESO_hatd(k); 
    RLS_regressor(k) = x(k);
    RLS_error(k) = RLS_output(k) - RLS_hata(k)*RLS_regressor(k);
    RLS_g(k) = RLS_P(k)*RLS_regressor(k)*( lambda + RLS_regressor(k)'*RLS_P(k)*RLS_regressor(k) )^(-1);
    RLS_P(k+1) = lambda^(-1)*RLS_P(k) - RLS_g(k)*RLS_regressor(k)'*lambda^(-1)*RLS_P(k);
    RLS_hata(k+1) = RLS_hata(k) + RLS_error(k)*RLS_g(k);
    
    %% ESO
    ESO_error(k) = x(k) - ESO_hatx(k);
    ESO_hatdx_dt(k) = RLS_hata(k)*x(k) + b_N*u_all(k) + ESO_hatd(k) + k1*ESO_error(k);
    ESO_hatdd_dt(k) = k2*ESO_error(k);
    ESO_hatx(k+1) = ESO_hatx(k) + dt*ESO_hatdx_dt(k);
    ESO_hatd(k+1) = ESO_hatd(k) + dt*ESO_hatdd_dt(k);
    
    %% Compute External Variables
    true_total_disturbance(k) = a(k)*x(k) + d(k) - RLS_hata(k)*x(k);
    hat_total_disturbance(k) = ESO_hatd(k); 
    
end

%% Prepare Plot
index_loop = 1:N+1; index_time = ( index_loop-1 )*dt; index_data = index_loop;

%% Plot Output
figure('name','System')
subplot(2,1,1); 
plot(index_time, xref(index_data), 'r-', 'linewidth', 1); hold on;
plot(index_time, x(index_data), 'b-.', 'linewidth', 1.5); grid on; 
title('System Output'); legend('x_{ref}', 'x')
subplot(2,1,2); 
plot(index_time, u_compensate(index_data), 'b-.', 'linewidth', 1.2); hold on;
plot(index_time, u_pid(index_data), 'r-.', 'linewidth', 1.2); hold on;
plot(index_time, u_all(index_data), 'k-', 'linewidth', 0.4); grid on; 
title('System Input'); legend('u_{compensate}', 'u_{pid}', 'u_{all}');
xlabel('Time [s]'); 

%% Plot RLS Results
figure('name','RLS')
subplot(2,1,1); 
plot(index_time, a(index_data), 'r-', 'linewidth', 1); hold on;
plot(index_time, RLS_hata(index_data), 'b-.', 'linewidth', 1.5); grid on;
title('Estimated Parameter'); legend('a', 'RLS_{hata}');
subplot(2,1,2); 
plot(index_time, RLS_error(index_data), 'r-', 'linewidth', 1); grid on;
title('Estimated Error'); 
xlabel('Time [s]');

%% Plot ESO Results
figure('name','ESO')
subplot(2,1,1); 
plot(index_time, x(index_data), 'r-', 'linewidth', 1); hold on;
plot(index_time, ESO_hatx(index_data),'b-.', 'linewidth', 1.5); grid on;
title('ESO Output'); legend('x', 'ESO_{hatx}');
subplot(2,1,2); 
plot(index_time, d(index_data), 'r-', 'linewidth', 1); hold on; 
plot(index_time, ESO_hatd(index_data), 'b-.', 'linewidth', 1.5); grid on; 
title('ESO Disturbance'); legend('d','ESO_{hatd}'); 
xlabel('Time [s]'); 

%% Plot Total Disturbance
figure('name','Total Disturbance')
plot(index_time, true_total_disturbance(index_data), 'r-', 'linewidth', 1); hold on;
plot(index_time, hat_total_disturbance(index_data), 'b-.', 'linewidth', 1.5); grid on;
title('Total Disturbance'); legend('true', 'hat');
xlabel('Time [s]');

%% Show the total run time
toc