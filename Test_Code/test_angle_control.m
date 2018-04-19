clc; clear; close all;
tic
% simulation parameters
dt = 0.0001; N = 600000;
% initial value
x(:,1) = [10,0]; y(1) = 5;
% system parameters
a = [3 4 5];
c = [2 3 4];
bx = 1; by = 5;
kpid_x = [3 0.000 5];
kpid_y = [1 0.000];
umax = 26;
for i = 1:N
    % reference
    yr(i) = -5*sin(dt*i);
    % controller
    error_x(i) = (-21*yr(i)/13)-x(1,i);
    error_y(i) = yr(i)-y(i);
    if i == 1
        error_x_sum(1) = (0-x(1,i))*dt;
        error_y_sum(1) = (0-y(i))*dt;
    else
        error_x_sum(i) = error_x_sum(i-1) + (0-x(1,i))*dt;
        error_y_sum(i) = error_y_sum(i-1) + (0-y(i))*dt;
    end
    u_x_pid(i) = kpid_x(1)*error_x(i) + kpid_x(2)*error_x_sum(i) + kpid_x(3)*(0-x(2,i));
    u_y_pid(i) = kpid_y(1)*error_y(i) + kpid_y(2)*error_y_sum(i);
    u(i) = u_x_pid(i) + u_y_pid(i);
    u(i) = min(max(u(i),-umax),umax);
    % system dynamics
    x(1,i+1) = x(1,i) + dt*( x(2,i) );
    x(2,i+1) = x(2,i) + dt*( a(1)*x(1,i)+a(2)*x(2,i)+a(3)*y(i)+bx*u(i) + 2*sin(dt*i) + 10*randn);
    y(i+1) = y(i) + dt*( c(1)*x(1,i)+c(2)*x(2,i)+c(3)*y(i)+by*u(i) + 2*cos(dt*i) + 10*randn);
end
% plot
index_time = dt*((1:N)-1);
index_data = 1:N;
figure('name','States')
subplot(2,2,1)
plot(index_time,u(index_data),'linewidth',1.2); legend('u'); xlabel('Time [s]'); grid on;
subplot(2,2,2)
plot(index_time,yr(index_data),'r-.','linewidth',0.8); hold on;
plot(index_time,y(index_data),'b-','linewidth',1.2); legend('yr','y'); xlabel('Time [s]'); grid on;
subplot(2,2,3)
plot(index_time,x(1,index_data),'linewidth',1.2); legend('x1'); xlabel('Time [s]'); grid on;
subplot(2,2,4)
plot(index_time,x(2,index_data),'linewidth',1.2); legend('x2'); xlabel('Time [s]'); grid on;
toc