clear all; clc; close all;

% Discretization of continuous system
dt = 0.01;
data = csvread('20_02_protocol/20_02_Protocol_Power_grasp_3.csv',1,0);
% Generate time array
data = data(2300:5804,:);
tpg = data(:,1)';

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
SignalNoisy = data(:,2:3)';

% LQE function
xp = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);figure
X = [SignalNoisy(1,:); SignalNoisy(2,:)]';
GMModel = fitgmdist(X,3);
figure
y = [zeros(1000,1);ones(1000,1)];
h = gscatter(X(:,1),X(:,2));
hold on
gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
axis([0 inf 0 5]);
hold off
%%
figure(1)
splitsPG = 60;
itpPG = findchangepts(xp(1,:),'MaxNumChanges',splitsPG,'Statistic','mean');
figure();set(gcf,'color','white');
hold on;
itptrim = [1];
for i= 1:splitsPG-1
    if itpPG(i+1)-itpPG(i) >60
        itptrim = [itptrim itpPG(i+1)];
    end
    
end
itpPG = itptrim;
plot(tpg,SignalNoisy(1,:),'r');
plot(tpg,SignalNoisy(2,:));
itpSize = size(itpPG);
for i=1:itpSize(2)
    xline(tpg(itpPG(i)));
end
hold off
%% Plots
set(0,'DefaultFigureWindowStyle','docked')

figure();set(gcf,'color','white');
hold on;

plot(tpg,SignalNoisy(1,:),'r');
plot(tpg,SignalNoisy(2,:));
plot(tpg,xp(1,:),'b','linewidth',1.5);
title('Kalman Filter 2 signals & sigma: 0.5')
legend('Sensor 1','Sensor 2','Combined Filtered Signal');

xlabel('Time (s)','fontsize',15);
ylabel('Voltage (V)','fontsize',15);



