clc; clear; close all
tic;
%% Simulation Setting
t_end = 20;
dt = 0.001;
Ts_sample = 0.01;
N_sample = floor(Ts_sample/dt);
t_log = 0:dt:t_end;
N = length(t_log);
% Ref
A_ref = 2;
T_ref = 2.5;
% STR
kp = 1;
% SMC
alpha = 2;
% RLS
forgettingfactor = 1;
hattheta(1) = 0.2;
P{1} = 1000;
% System
dp_exp = 1.0;
b = 1.2;
bN = 1.0;
x(1) = 5;
% Key Log
k = 0;
t_sample = [];
state_switch = [];
%% Start Simulation
for i = 1:N
    t_now = dt*(i-1);
    theta(i) = 2 + 0.2*sin(2*pi/1*t_now);
    dp(i) = 0.01 + t_now^(dp_exp);
    %% Sampling  and Computing
    if ( mod(i,N_sample) == 1 )
        k = k + 1;
        %% Sampling
        t_sample = [t_sample,t_now];
        theta_true(k) = theta(i);
        dp_true(k) = dp(i);
        dp_measure(k) = dp_true(k) + 0.00*randn;
        if i<=10
            x_true(k) = x(i);
        else
            x_true(k) = x(i-10);
        end
        x_measure(k) = x_true(k) + 0.000*randn;
        %% RLS
        if k == 1
            error_RLS(k) = 0;
            temp = P{k};
            P_diag(1,k) = temp(1,1);
        else
            output(k) = ( x_measure(k) - x_measure(k-1) ) / Ts_sample - bN*u_design( k-1 );
            regressor(k-1) = x_measure(k-1)*dp_measure(k-1);
            error_RLS(k) = output(k) - hattheta(k-1)*regressor(k-1);
            K(k) = P{k-1}*regressor(k-1)/(forgettingfactor + regressor(k-1)'*P{k-1}*regressor(k-1));
            if (abs(error_RLS(k)))<2  % error based 
                P{k} = forgettingfactor^(-1)*P{k-1}-K(k)*regressor(k-1)'*forgettingfactor^(-1)*P{k-1};
            else
                P{k} = 1000;
            end
            temp = P{k};
            P_diag(1,k) = temp(1,1);
            hattheta(k) = hattheta(k-1) + max(min(K(k)*error_RLS(k),1e-4+0.02*abs(hattheta(k-1))),-1e-4-0.02*abs(hattheta(k-1)));
        end
        %% ESO
        pole = 15;
        if k == 1 
            hatx(k) = x_measure(k);
            hatd(k) = 0;
        else
            hatx(k) = hatx(k-1) + Ts_sample*( hattheta(k)*x_measure(k-1)*dp_measure(k) + bN*u_design(k-1) + hatd(k-1) - 2*pole*(hatx(k-1) - x_measure(k)) );
            hatd(k) = hatd(k-1) + Ts_sample*( -pole^2*(hatx(k-1)-x_measure(k)) );
        end
        %% Reference
        x_ref(k) = A_ref*sin(2*pi/T_ref*t_now);
        dx_ref(k) = 2*pi/T_ref*A_ref*cos(2*pi/T_ref*t_now);
        %% Controller
        error_tracking(k) = x_measure(k) - x_ref(k);
        if abs(error_tracking(k))<1
            u_design(k) = ( dx_ref(k) - 5*error_tracking(k) - hattheta(k)*x_measure(k)*dp_measure(k) - 0*hatd(k))/bN;
            state_switch = [state_switch,0];
            if k == 1
                Gamma(1) = 50;
            else
                Gamma(k) = Gamma(k-1);
            end
        else
            sign_error_tracking = ( 1-exp(-alpha*error_tracking(k)) ) / ( 1+exp(-alpha*error_tracking(k)) );
            if k == 1
                Gamma(1) = 50;
            else
                Gamma(k) = Gamma(k-1) + Ts_sample*0.1*error_tracking(k)*sign_error_tracking;
            end
            u_design(k) = ( dx_ref(k) - Gamma(k)*sign_error_tracking - hattheta(k)*x_measure(k)*dp_measure(k) - 0*hatd(k))/bN;
            state_switch = [state_switch,1];
        end
        u_design(k) = ( dx_ref(k) - 5*error_tracking(k) - hattheta(k)*x_measure(k)*dp_measure(k) - 1*hatd(k) )/bN;
    end
    %% System
    x(i+1) = x(i) + dt*( theta(i)*x(i)*dp(i) + b*u_design(k) );
end
%% Plot
index_plot = 1:length(t_sample);
figure('name','Estimator')
subplot(3,1,1)
plot(t_sample,error_RLS(index_plot),'linewidth',2); hold on;
grid on;
legend({'e'})
subplot(3,1,2)
plot(t_sample,P_diag(1,index_plot),'r-','linewidth',2); hold on;
grid on;
legend({'P11'})
subplot(3,1,3)
plot(t_sample,theta(index_plot),'r-','linewidth',3); hold on;
plot(t_sample,hattheta(index_plot),'b-','linewidth',1.5); hold on;
grid on;
legend('True','RLS');
ylim([-2,6]);
xlabel('Time [s]')
figure('name','TrackingResults')
plot(t_sample,x_ref(index_plot),'r-','linewidth',3); hold on;
plot(t_sample,x_true(index_plot),'b-.','linewidth',1.5); hold on;
grid on;
legend({'x_{ref}','x'});
xlabel('Time [s]')
figure('name','switch state')
plot(t_sample,state_switch(index_plot),'bo','linewidth',2); hold on;
xlabel('Time [s]');
toc