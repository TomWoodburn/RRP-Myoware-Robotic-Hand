%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LQE - runs a Kalman Filter for one time loop.
% Inputs:
% A - State transition matrix
% C - Observation matrix
% Q - System noise covariance matrix
% R - Observation noise covariance matrix
% z0 - Initial value of state
% P0 - Initial value of error covariance matrix
% y - Noisy observations
%
% Output:
% z - State variable over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [z] = KalmanFilter(A,C,Q,R,z0,P0,y)

    % Allocate memory for speed
    z = zeros(size(z0,1),length(y));
    % Initialize Kalman filter
    z(:,1) = z0;
    P = P0;

    % Kalman filter loop
    for i=1:length(y)-1
        % Prediction
        z_Prior = A*z(:,i);
        P = A*P*A'+Q;

        % Correction
        K = P*C'/(C*P*C'+R);
        z(:,i+1) = z_Prior + K*(y(:,i+1)-C*z_Prior);
        P = (eye(3)-K*C)*P;
    end