clc; clear;close all;
tic;
% Define the data to be fit
alpha = (-25:5:25)';        % angle of attack
beta = (-25:5:25)';         % angle of slide
delta_phi = (-30:10:30)';   % deflection of pitch dynamics
delta_psi = (-30:10:30)';   % deflection of yaw dynamics
delta_gamma = (-30:10:30)'; % deflection of roll dynamics
CN_zero = zeros(length(alpha),length(beta));    % part I
CN_pitch = zeros(length(alpha),length(beta),length(delta_phi)); % part II
CN_yaw = zeros(length(alpha),length(beta),length(delta_psi));   % part III
CN_roll = zeros(length(alpha),length(beta),length(delta_gamma));% part IV
scale = 0.1;   % noise

%%
%   Cmx = (-p1/Ma^2+p2+p3*alpha+p4*beta)*alpha + ...
%         (-p5/Ma^2+p6+p4*alpha-p3*beta)*beta+ ...
%         (p7*int_{0,t}dot{w_{xb}}+p8*int_{0,t}dot{w_{yb}}+p9*int_{0,t}dot{w_{zb}})
%
for i_alpha = 1:length(alpha)
    for i_beta = 1:length(beta)
        CN_zero(i_alpha,i_beta) = 2e-8 + 1.5*alpha(i_alpha) + 2e-6*beta(i_beta) + 0.5*cos(alpha(i_alpha)*pi/180) + scale*randn;
        for i_delta_phi = 1:length(delta_phi)
            CN_pitch(i_alpha,i_beta,i_delta_phi) = 2e-9 + 0.5*alpha(i_alpha) + 2e-7*beta(i_beta) - 0.6*delta_phi(i_delta_phi) + scale*randn;
        end
        for i_delta_psi = 1:length(delta_psi)
            CN_yaw(i_alpha,i_beta,i_delta_psi) = 2e-10 + 0.2*alpha(i_alpha) + 2e-5*beta(i_beta) + 2e-4*delta_psi(i_delta_psi) + scale*randn;
        end
        for i_delta_gamma = 1:length(delta_gamma)
            CN_roll(i_alpha,i_beta,i_delta_gamma) = 2e-9 + 0.3*alpha(i_alpha) + 2e-8*beta(i_beta) + 2e-5*delta_gamma(i_delta_gamma) + scale*randn;
        end
    end
end
alpha_fit = zeros(11*11*7*7*7,1);
delta_phi_fit = zeros(11*11*7*7*7,1);
CN = zeros(11*11*7*7*7,1);
for i_alpha = 1:length(alpha)
    for i_beta = 1:length(beta)
        for i_delta_phi = 1:length(delta_phi)
            for i_delta_psi = 1:length(delta_psi)
                for i_delta_gamma = 1:length(delta_gamma)
                    CN_zero_chosen = CN_zero(alpha==alpha(i_alpha),beta==beta(i_beta));
                    CN_pitch_chosen = CN_pitch(alpha==alpha(i_alpha),beta==beta(i_beta),delta_phi==delta_phi(i_delta_phi));
                    CN_yaw_chosen = CN_yaw(alpha==alpha(i_alpha),beta==beta(i_beta),delta_psi==delta_psi(i_delta_psi));
                    CN_roll_chosen = CN_roll(alpha==alpha(i_alpha),beta==beta(i_beta),delta_gamma==delta_gamma(i_delta_gamma));
                    alpha_fit(11*7*7*7*(i_alpha-1)+7*7*7*(i_beta-1)+7*7*(i_delta_phi-1)+7*(i_delta_psi-1)+i_delta_gamma)...
                        =alpha(i_alpha);
                    delta_phi_fit(11*7*7*7*(i_alpha-1)+7*7*7*(i_beta-1)+7*7*(i_delta_phi-1)+7*(i_delta_psi-1)+i_delta_gamma)...
                        =delta_phi(i_delta_phi);
                    CN(11*7*7*7*(i_alpha-1)+7*7*7*(i_beta-1)+7*7*(i_delta_phi-1)+7*(i_delta_psi-1)+i_delta_gamma)...
                        =CN_zero_chosen + CN_pitch_chosen + CN_yaw_chosen + CN_roll_chosen;
                end
            end
        end
    end
end
% % Define function that will be used to fit data
% % (P is a vector of fitting parameters)
f = @(P,x) P(1) + P(2)*x(:,1) + P(3)*x(:,2) + P(4)*cos(x(:,1)*pi/180);
[P_fitted,R,J,CovB,MSE,ErrorModelInfo] = nlinfit([alpha_fit,delta_phi_fit],CN,f,[1 1 1,1]);
% % Display fitted coefficients
disp(['P = ',num2str(P_fitted)])
% % Plot the data and fit
figure
plot(CN,'*'); hold on;
plot(f(P_fitted,[alpha_fit,delta_phi_fit]),'g');
legend('data','fit')
toc