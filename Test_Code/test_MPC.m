%%  Clear
clear; close all;
tic

%%  Simulation setting
dt = 0.001; % 1ms
Ts = 0.01;  % 10ms
t_end = 60;
t_log = 0:dt:t_end;
t_idx_all = length(t_log);

%%  Construct model
A_c = [-10 0;
    90 -10];
B_c = [10;
    20];
A_d = eye(2)+Ts*A_c;
B_d = Ts*B_c;
C_d = [1 1];
D_d = 0;
n_in = size(B_d,2);
n_states = size(A_d,1);
n_out = size(C_d,1);

N_mpc = 50;
Q_mpc = 1*eye(n_out);
S_mpc = 1000*eye(n_in);



%%	Initial matrix filling
Lambda_mpc = zeros(N_mpc*n_out,n_states);
for i = 1:N_mpc
    Lambda_mpc((i-1)*n_out+(1:n_out),:) = C_d*A_d^(i-1);
end
AB_mpc = zeros(N_mpc*n_out,n_in);
AB_mpc(1:n_out,:) = D_d;
for i = 2:N_mpc
    AB_mpc((i-1)*n_out+(1:n_out),:) = C_d*A_d^(i-2)*B_d;
end
Phi_mpc = zeros(N_mpc*n_out);
for i = 1:N_mpc
    Phi_mpc(1+(i-1)*n_out:end,(i-1)*n_in+(1:n_in)) = AB_mpc(1:(N_mpc-i+1)*n_out,:);
end
barQ_mpc = kron(eye(N_mpc),Q_mpc);
barS_mpc = kron(diag(2*ones(1,N_mpc),0),S_mpc) ...
    + kron(diag(-ones(1,N_mpc-1), 1),S_mpc) ...
    + kron(diag(-ones(1,N_mpc-1),-1),S_mpc);
H_mpc = Phi_mpc'*barQ_mpc*Phi_mpc + barS_mpc;

%%	Create an observer
Aa = blkdiag(A_d,eye(n_out)); Ba = [B_d;zeros(n_out,n_in)];
Ca = [C_d eye(n_out)]; Da = D_d;
Ga = diag([ones(1,n_states),ones(1,n_out)]);
Qa = 1*eye(n_states+n_out);
Ra = 1*eye(n_out);
L = dlqe(Aa,Ga,Ca,Qa,Ra);
L = Aa*L;

%%	Create plant constraints
lb_mpc = -0.1*ones(n_in,1); 
ub_mpc = 0.1*ones(n_in,1);
Lb_mpc = kron(ones(N_mpc,1),lb_mpc); 
Ub_mpc = kron(ones(N_mpc,1),ub_mpc);

%%  Preallocation memory
ref = zeros(n_out,t_idx_all);
x = zeros(n_states,t_idx_all);
hat_x = zeros(n_states+n_out,t_idx_all);   % Estimated state
u = zeros(n_in,t_idx_all);              % Input record
y = zeros(n_out,t_idx_all);             % Output record

%%  Set system initial conditions
x(:,1) = 0*ones(n_states,1);   % Plant state
u_old = zeros(n_in,1);
Aobs = Aa-L*Ca;     % Generate observer state gain matrix
Bobs = Ba-L*Da;     % Generate observer input gain matrix
fd_mpc = [Phi_mpc'*barQ_mpc*Lambda_mpc, Phi_mpc'*barQ_mpc];
opt = optimset('display','off');

%%  Generate noise signal
noise = 1e-3*cumsum(randn(t_idx_all,n_out)')' + 1e-3*randn(t_idx_all,n_out);
[bfilt,afilt] = butter(4,0.9);
noise = filter(bfilt,afilt,noise)';

%% Start simulation
for i = 1:t_idx_all
    
    %% time now
    t_now = dt*(i-1);
    
    %% true system
    x(:,i+1) = x(:,i) + dt*( A_c*x(:,i) + B_c*u(:,i)+ 0.2*sin(t_now)*ones(n_states,1) );
    y(:,i) = C_d*x(:,i) + D_d*u(:,i) + noise(:,i);
    
    %% sampling
    if mod(t_now,Ts) == 0
        %% reference
        if t_now<=t_end/4
            ref(:,i) = 1;
        elseif t_now<=t_end*2/4
            ref(:,i) = -1;
        elseif t_now<=t_end*3/4
            ref(:,i) = 1;
        else
            ref(:,i) = -1;
        end
        %% controller
        hat_x(:,i+floor(Ts/dt)) = Aobs*hat_x(:,i) + Bobs*u(:,i) + L*y(:,i);
        ref_mpc = kron(ones(N_mpc,1),ref(:,i));
        f = Phi_mpc'*barQ_mpc*( Lambda_mpc*hat_x(1:n_states,i+floor(Ts/dt))-ref_mpc ) ...
            - [S_mpc*u_old;zeros((N_mpc-1)*n_out,1)];
        U = quadprog(H_mpc,f,[],[],[],[],Lb_mpc,Ub_mpc,[],opt);
        u(:,i+floor(Ts/dt)) = U(1:n_in); u_old = u(:,i);
    else
        ref(:,i) = ref(:,i-1);
        hat_x(:,i+floor(Ts/dt)) = hat_x(:,i-1);
        u(:,i+floor(Ts/dt)) = u(:,i-1);
    end
end

%% Plot
time_plot = t_log;
idx_plot = 1:t_idx_all;
figure('name','Control Results')
subplot(2,1,1);
plot(time_plot,ref(idx_plot)','r--','Linewidth',4); hold on
plot(time_plot,(x(1,idx_plot)+x(2,idx_plot))','b-','Linewidth',2); hold on;
plot(time_plot,y(idx_plot)','m-','Linewidth',2); hold on;
grid on;
legend({'Ref','x','y'});
subplot(2,1,2);
plot(time_plot,ub_mpc*ones(1,t_idx_all),'r--','Linewidth',4); hold on;
plot(time_plot,lb_mpc*ones(1,t_idx_all),'k--','Linewidth',4); hold on;
plot(time_plot,u(:,idx_plot)','b-','Linewidth',2); hold off;
grid on;
xlabel('Time [s]');
legend({'u_{max}','u_{min}','u'});

toc