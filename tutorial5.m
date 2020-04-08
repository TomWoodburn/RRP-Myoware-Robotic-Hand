clear all; clc; close all;

% Discretization of continuous system
dt = 0.01;
data = csvread('20_02_protocol\20_02_Protocol_Power_grasp_1_s.csv',1,0);
% Generate time array

t = data(:,1)';

%%%%%%% Kalman filter matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State transition matrix
A = [1 dt dt^2/2;
     0  1     dt;
     0  0     1];

% System noise and covariance matrix
SigmaQ = 0.3;
Q=[ SigmaQ^6/36    SigmaQ^5/12   SigmaQ^4/6
    SigmaQ^5/12    SigmaQ^4/4    SigmaQ^3/2
    SigmaQ^4/6     SigmaQ^3/2    SigmaQ^2];

% Observation matrix
C = [1 0 0;
    1 0 0];

% Observation noise and covariance matrix
Sigma = 0.5;
R=[(Sigma^2) 0;
    0 (Sigma^2)];

% Initialize state and error covariance matrices
xInit = zeros(3,1);
PInit = diag([1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate signal and noise corrupted signal
SignalNoisy = data(:,2:3)';

% LQE function
 x = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);

%% Plots
set(0,'DefaultFigureWindowStyle','docked')

figure(1);set(gcf,'color','white');
hold on;

plot(t,SignalNoisy(1,:),'r');
plot(t,SignalNoisy(2,:));
plot(t,x(1,:),'b','linewidth',1.5);
title('Kalman Filter 2 signals & sigma: 0.5')
legend('Sensor 1','Sensor 2','Combined Filtered Signal');

xlabel('Time (s)','fontsize',15);
ylabel('Voltage (V)','fontsize',15);



