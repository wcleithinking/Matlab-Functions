%% test_EKF.m
%--------------------------------------------------------------------------
%   System: 
%       \dox{x} = f_sys(x) + noise1
%       y = h_sys(x) + noise2
%--------------------------------------------------------------------------
clc; clear; close all;
tic;
%% Simulation Setting
t_end = 60;
dt = 0.02;
t_log = 0:dt:t_end;
N = length(t_log);
%% Initial Values
x(1) = 10;
y(1) = x(1) + 2*randn;
hatx(1) = 0;
P{1} = 100;
%% Start Simulation
for k = 1:N
    t_now = dt*(k-1);
    % Discrete-Time System
    x(k+1) = f_sys(x(k)) + 0.2*randn;
    y(k+1) = h_sys(x(k+1)) + 2*randn;
    % EKF
    [hatx(k+1), P{k+1}] = myEKF(@f_sys, @df_sys, @h_sys, @dh_sys, hatx(k), P{k}, y(k), 0.2^2, 2^2);
end
%% Plot
index_plot = 1:N;
plot(t_log,y(index_plot),'b-.','linewidth',0.5); hold on;
plot(t_log,hatx(index_plot),'g-.','linewidth',2); hold on;
plot(t_log,x(index_plot),'r-','linewidth', 3); hold on;
grid on;
legend('y','EKF','True')
xlabel('Time [s]')
toc