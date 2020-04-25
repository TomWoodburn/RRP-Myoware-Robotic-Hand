clear all; clc; close all;
set(0,'DefaultFigureWindowStyle','docked')
% Discretization of continuous system
dt = 0.01;
dataS = csvread('20_02_protocol/20_02_Protocol_sphere_finger_1.csv',1,0);
% Generate time array
dataS = dataS(800:4700,:);
ts = dataS(:,1)';

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
R=[(0.5^2) 0;
    0 (0.5^2)];

% Initialize state and error covariance matrices
xInit = zeros(3,1);
PInit = diag([1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate signal and noise corrupted signal
SignalNoisy = dataS(:,2:3)';

% LQE function
xs = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);
figure
X = [SignalNoisy(1,:); SignalNoisy(2,:)]'
GMModel = fitgmdist(X,2);
figure
y = [zeros(1000,1);ones(1000,1)];
h = gscatter(X(:,1),X(:,2));
hold on
gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
axis([0 inf 0 5])
hold off
%%
figure(1)
splitsS = 60;
itpS = findchangepts(xs(1,:),'MaxNumChanges',splitsS,'Statistic','mean');
figure();set(gcf,'color','white');
hold on;
itptrim = [1];
for i= 1:splitsS-1
    if itpS(i+1)-itpS(i) >60
        itptrim = [itptrim itpS(i+1)];
    end
    
end
itpS = itptrim;
plot(ts,SignalNoisy(1,:),'r');
plot(ts,SignalNoisy(2,:));
itpSize = size(itpS);
for i=1:itpSize(2)
    xline(ts(itpS(i)));
end
hold off

%% Plots


figure();set(gcf,'color','white');
hold on;

plot(ts,SignalNoisy(1,:),'r');
plot(ts,SignalNoisy(2,:));
plot(ts,xs(1,:),'b','linewidth',1.5);
title('Kalman Filter 2 signals & sigma: 0.5')
legend('Sensor 1','Sensor 2','Combined Filtered Signal');

xlabel('Time (s)','fontsize',15);
ylabel('Voltage (V)','fontsize',15);



