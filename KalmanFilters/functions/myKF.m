function [hatx_new, P_new] = myKF(A_sys, C_sys, hatx_old, P_old, y_new, Q, R)
%--------------------------------------------------------------------------
% This function provides the single-step iterative process of the standard Kalman Filter (KF).
% Author: Wenchao Lei.
% Date:   2018-04-15
%   DESCRIPTION:
%       Inputs:
%           A_sys:    the state transition matrix of the discrete-time system
%           C_sys:    the measurement matrix of the discrete-time system
%           hatx_old: the old estimate of state x
%           P_old:    the old covariance matrix under estimate hatx_old
%           y_new:    the new measurement of output y 
%           Q:        the covariance matrix of the process noise
%           R:        the covariance matrix of the measurement noise
%       Outputs:
%           hatx_new: the new estimate of state x
%           P_new:    the new covariance matrix under estimate hatx_new
%--------------------------------------------------------------------------

%-------------------------------- Main -----------------------------------
% predict
prehatx_new = A_sys*hatx_old;
preP_new = A_sys*P_old*A_sys' + Q;
% update
preerror = y_new - C_sys*prehatx_new;
S = C_sys*preP_new*C_sys' + R;
K = preP_new*C_sys/S;
hatx_new = prehatx_new + K*preerror;
P_new = ( eye(length(hatx_old)) - K*C_sys )*preP_new;
end