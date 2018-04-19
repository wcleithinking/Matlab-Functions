clc
close
clear
%% Preprare Simulation
% simulation parameters
dt = 0.0005;
N = 120000;
% system parameters
b = 2;
% initial data
x(1) = 10;
output(1) = 0;
regressor(1) = 0;
hata(1) = 0; 
hatx(1) = 0;
hatdis(1) = 0;
error_x(1) = 0;
% tunable parameters
lambda = 0.985;
P(1) = 100;
w_tuning = 5;
k1 = 2*w_tuning;
k2 = w_tuning^2;
for k = 1:N+1
    %% Controller
    u(k) = min(max((-x(k)-hata(k)*x(k)-hatdis(k))/b,-16),16);
    %% System
    a(k) = 2.5+0.5*rand+0.5*sin(dt*k);
    dis(k) = 1.5*sin(dt*k)+1;
    dx(k) = a(k)*x(k) + b*u(k) + dis(k);
    x(k+1) = x(k) + dt*dx(k);
    %% RLS
    error_RLS(k) = output(k) - hata(k)*regressor(k);
    g(k) = P(k)*regressor(k)*(lambda + regressor(k)'*P(k)*regressor(k))^(-1);
    P(k+1) = lambda^(-1)*P(k)-g(k)*regressor(k)'*lambda^(-1)*P(k);
    hata(k+1) = hata(k) + error_RLS(k)*g(k);
    output(k+1) = ( x(k+1)-x(k) )/dt - b*u(k) - hatdis(k);
    regressor(k+1) = x(k);
    %% ESO
    hatdx(k) = hata(k)*x(k) + b*u(k) + hatdis(k) + k1*error_x(k);
    hatddis(k) = k2*error_x(k);
    hatx(k+1) = hatx(k) + dt*hatdx(k);
    hatdis(k+1) = hatdis(k) + dt*hatddis(k);
    error_x(k+1) = x(k+1) - hatx(k+1);
end
%% Prepare Plot
index_time = ((1:N+1)-1)*dt;
index_data = 1:N+1;
%% Plot Output
figure('name','x')
subplot(2,1,1)
plot(index_time,x(index_data),'b-'); grid on;
title('x')
subplot(2,1,2)
plot(index_time,u(index_data),'k-'); grid on;
title('u');
xlabel('Time [s]'); 
%% Plot RLS
figure('name','RLS')
subplot(2,1,1)
plot(index_time, a(index_data),'r-',index_time, hata(index_data),'b-.');
title('Estimated Parameter'); legend('a','hata');
subplot(2,1,2)
plot(index_time, error_RLS(index_data),'r-');
title('Error')
xlabel('Time [s]');
%% Plot ESO
figure('name','ESO')
subplot(2,1,1)
plot(index_time, x(index_data),'r-',index_time, hatx(index_data),'b-.'); grid on;
ylabel('x');legend('x','hatx');
subplot(2,1,2)
plot(index_time, dis(index_data),'r-',index_time, hatdis(index_data),'b-.'); grid on;
ylabel('Disturbance');legend('dis','hatdis');
xlabel('Time [s]'); 