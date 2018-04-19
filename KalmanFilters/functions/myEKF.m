function [hatx_new, P_new] = myEKF(f_sys, df_sys, h_sys, dh_sys, hatx_old, P_old, y_new, Q, R)
%--------------------------------------------------------------------------
% This function provides the single-step iterative process of the Extended Kalman Filter (EKF).
% Author: Wenchao Lei.
% Date:   2018-04-15
%   DESCRIPTION:
%       Inputs:
%           f_sys:    the f function of the discrete-time system
%           h_sys:    the h function of the discrete-time system
%           df_sys:   the derivative function of f function
%           dh_sys:   the derivative function of h function
%           hatx_old: the old estimate of state x
%           P_old:    the old covariance matrix under estimate hatx_old
%           y_new:    the new measurement of output y 
%           Q:        the covariance matrix of the process noise
%           R:        the covariance matrix of the measurement noise
%       Outputs:
%           hatx_new: the new estimate of state x
%           P_new:    the new covariance matrix under estimate hatx_new
%--------------------------------------------------------------------------

%-------------------------------- Main ------------------------------------
% predict
F = df_sys(hatx_old);
prehatx_new = f_sys(hatx_old);
preP_new = F*P_old*F' + Q;
% update
H = dh_sys( prehatx_new );
preerror = y_new - h_sys( prehatx_new );
S = H*preP_new*H' + R;
K = preP_new*H/S;
hatx_new = prehatx_new + K*preerror;
P_new = ( eye(length(hatx_old)) - K*H )*preP_new;