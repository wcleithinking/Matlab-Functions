clc; clear; close all; format long;
tic
dt = 0.001;
N = 6000;
x1(1) = 25;
y1(1) = 25;
x2(1) = 25;
y2(1) = 25;
for t=1:N
    % xref system
    x1ref(t) = 20*sin(dt*t*5);
    x2ref(t) = 20*cos(dt*t*5);
    y1ref(t) = 20*sin(dt*t*5);
    y2ref(t) = 20*cos(dt*t*5);
    % x system
    ex1(t) = x1ref(t) - x1(t);
    ex2(t) = x2ref(t) - x2(t);
    kx1(t) = 5;
    kx2(t) = 5;
    ux(t) = kx1(t)*ex1(t)+ kx2(t)*ex2(t);
    x1(t+1) = x1(t) + dt*( 10*x2(t)+5*sin(dt*t*2) + 2*ux(t) );
    x2(t+1) = x2(t) + dt*( 2*x1(t)+5*cos(dt*t*2) + 5*ux(t) );
    % y system
    ey1(t) = y1ref(t) - y1(t);
    ey2(t) = y2ref(t) - y2(t);
    ky1(t) = 5*(abs(ey1(t)))^2 / ( (abs(ey1(t)))^2+(abs(ey2(t)))^2 );
    ky2(t) = 5*(abs(ey2(t)))^2 / ( (abs(ey1(t)))^2+(abs(ey2(t)))^2 );
    uy(t) = ky1(t)*ey1(t)+ ky2(t)*ey2(t);
    y1(t+1) = y1(t) + dt*( 10*y2(t)+5*sin(dt*t*2) + 2*uy(t) );
    y2(t+1) = y2(t) + dt*( 2*y1(t)+5*cos(dt*t*2) + 5*uy(t) );
end
time_index = ( (1:N) - 1 )*dt;
data_index = 1:N;
figure(1)
subplot(3,1,1)
plot(time_index, zeros(N,1), 'r-.'); hold on;
plot(time_index, ux(data_index), 'b-'); hold on;
plot(time_index, uy(data_index), 'k-'); hold on;
grid on;
subplot(3,1,2)
plot(time_index, x1ref(data_index), 'r-.'); hold on;
plot(time_index, x1(data_index), 'b-'); hold on;
plot(time_index, y1(data_index), 'k-'); hold on;
grid on;
subplot(3,1,3)
plot(time_index, x2ref(data_index), 'r-.'); hold on;
plot(time_index, x2(data_index), 'b-'); hold on;
plot(time_index, y2(data_index), 'k-'); hold on;
grid on;
toc