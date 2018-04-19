%% test_KF.m
%--------------------------------------------------------------------------
%   System: 
%       \dox{x} = -0.2x + noise1
%       y = x + noise2
%--------------------------------------------------------------------------
clc; clear; close all;
tic;
%% Simulation Setting
t_end = 60;
dt = 0.02;
t_log = 0:dt:t_end;
N = length(t_log);
A = 1-0.2*dt;
C = 1;
%% Initial Values
x(1) = 10;
y(1) = x(1) + 2*randn;
hatx(1) = 0;
P{1} = 100;
%% Start Simulation
for k = 1:N
    t_now = dt*(k-1);
    % Discrete-Time System
    x(k+1) = A*x(k) + 0.2*randn;
    y(k+1) = C*x(k+1) + 2*randn;
    % KF
    [hatx(k+1), P{k+1}] = myKF(A, C, hatx(k), P{k}, y(k), 0.2^2, 2^2);
end
%% Plot
index_plot = 1:N;
plot(t_log,y(index_plot),'b-.','linewidth',0.5); hold on;
plot(t_log,hatx(index_plot),'g-.','linewidth',2); hold on;
plot(t_log,x(index_plot),'r-','linewidth', 3); hold on;
grid on;
legend('y','KF','True')
xlabel('Time [s]')
toc