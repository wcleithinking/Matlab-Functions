function [hatx_new, P_new] = myUKF(f_sys, h_sys, hatx_old, P_old, y_new, E_nx, E_ny, Q, R, alpha, beta, kappa)
%--------------------------------------------------------------------------
% This function provides the single-step iterative process of the Unscented Kalman Filter (UKF).
% For more details, please see https://en.wikipedia.org/wiki/Kalman_filter#Unscented_Kalman_filter
% Author: Wenchao Lei.
% Date:   2018-04-16
%   DESCRIPTION:
%       Inputs:
%           f_sys:    the f function in the discrete-time system
%           h_sys:    the h function in the discrete-time system
%           hatx_old: the old estimate of state x
%           P_old:    the old covariance matrix under estimate hatx_old
%           y_new:    the new measurement of output y 
%           E_nx:     the expectation of the process noise
%           E_ny:     the expectation of the measurement noise
%           Q:        the covariance matrix of the process noise
%           R:        the covariance matrix of the measurement noise
%           alpha:    the parameter control the spread of the sigma point,
%           beta:     the parameter related to the distribution of state
%           kappa:    the parameter control the spread of the sigma point
%       Outputs:
%           hatx_new: the new estimate of state x
%           P_new:    the new covariance matrix under estimate hatx_new
%--------------------------------------------------------------------------

%-------------------------------- Main ------------------------------------
% predict

hatx_old_aug = [hatx_old; E_nx];
P_old_aug = blkdiag(P_old,Q)+1e-6;

L_pre = length(hatx_old_aug);
lambda_pre = alpha^2*( L_pre + kappa ) - L_pre;

sqrtP_old_aug = chol(P_old_aug,'lower');

prehatx_new = lambda_pre/(L_pre+lambda_pre)*f_sys( hatx_old_aug(1:length(hatx_old)) );
for i=1:2*L_pre
    
    if i<=L_pre
        sp_pre(1:L_pre,i) = hatx_old_aug + sqrt(L_pre+lambda_pre)*sqrtP_old_aug(1:L_pre,i);
    else
        sp_pre(1:L_pre,i) = hatx_old_aug - sqrt(L_pre+lambda_pre)*sqrtP_old_aug(1:L_pre,i-L_pre);
    end
    
    sp_process(1:length(hatx_old),i) = f_sys( sp_pre(1:length(hatx_old),i) );
    prehatx_new = prehatx_new + 1/(2*(L_pre+lambda_pre))*sp_process(1:length(hatx_old),i);
    
end

preP_new = ( lambda_pre/(L_pre+lambda_pre)+(1-alpha^2+beta) )*( f_sys( hatx_old_aug(1:length(hatx_old)) ) - prehatx_new )*( f_sys( hatx_old_aug(1:length(hatx_old)) ) - prehatx_new )';
for i=1:2*L_pre
    preP_new = preP_new + 1/(2*(L_pre+lambda_pre))*( sp_process(1:length(hatx_old),i) - prehatx_new )*( sp_process(1:length(hatx_old),i) - prehatx_new )';
end

% update

prehatx_new_aug = [prehatx_new; E_ny];
preP_new_aug = blkdiag(preP_new,R)+1e-8;

L_upd = length(prehatx_new_aug);
lambda_upd = alpha^2*( L_upd + kappa ) - L_upd;

sqrtpreP_new_aug = chol(preP_new_aug,'lower');

sp_upd0 = prehatx_new_aug;
sp_measurement0 = h_sys( sp_upd0(1:length(prehatx_new)) );
haty_new = lambda_upd/(L_upd+lambda_upd)*sp_measurement0;
for i=1:2*L_upd
    
    if i<=L_upd
        sp_upd(1:L_upd,i) = prehatx_new_aug + sqrt(L_upd+lambda_upd)*sqrtpreP_new_aug(1:L_upd,i);
    else
        sp_upd(1:L_upd,i) = prehatx_new_aug - sqrt(L_upd+lambda_upd)*sqrtpreP_new_aug(1:L_upd,i-L_upd);
    end
    
    sp_measurement(1:length(y_new),i) = h_sys( sp_upd(1:length(prehatx_new),i) );
    haty_new = haty_new + 1/(2*(L_upd+lambda_upd))*sp_measurement(1:length(y_new),i);
    
end

Pyy_new = ( lambda_upd/(L_upd+lambda_upd)+(1-alpha^2+beta) )*( sp_measurement0 - haty_new )*( sp_measurement0 - haty_new )';
Pxy_new = ( lambda_upd/(L_upd+lambda_upd)+(1-alpha^2+beta) )*( sp_upd0(1:length(prehatx_new)) - prehatx_new )*( sp_measurement0 - haty_new )';
for i=1:2*L_upd
    Pyy_new = Pyy_new + 1/(2*(L_upd+lambda_upd))*( sp_measurement(1:length(y_new),i) - haty_new  )*( sp_measurement(1:length(y_new),i) - haty_new )';
    Pxy_new = Pxy_new + 1/(2*(L_upd+lambda_upd))*( sp_upd(1:length(hatx_old),i) - prehatx_new )*( sp_measurement(1:length(y_new),i) - haty_new )';
end

K = Pxy_new/Pyy_new;
hatx_new = hatx_old + K*(y_new-haty_new);
P_new = preP_new - K*Pxy_new';