clear all; clc; close all;

% Discretization of continuous system
dt = 0.01;
data = csvread('20_02_protocol/20_02_Protocol_parallel_extension_1.csv',1,0);
% Generate time array
data = data(1050:4400,:);
tpe = data(:,1)';

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
xpe = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);
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
splitsPE = 90;
itpPE = findchangepts(xpe(1,:),'MaxNumChanges',splitsPE,'Statistic','mean');
figure();set(gcf,'color','white');
hold on;
itptrim = [1];
for i= 1:splitsPE-1
    if itpPE(i+1)-itpPE(i) >50
        itptrim = [itptrim itpPE(i+1)];
    end
    
end
itpPE = itptrim;
plot(tpe,SignalNoisy(1,:),'r');
plot(tpe,SignalNoisy(2,:));
itpSize = size(itpPE);
for i=1:itpSize(2)
    xline(tpe(itpPE(i)));
end
hold off
%% Plots
set(0,'DefaultFigureWindowStyle','docked')

figure();set(gcf,'color','white');
hold on;

plot(tpe,SignalNoisy(1,:),'r');
plot(tpe,SignalNoisy(2,:));
plot(tpe,xpe(1,:),'b','linewidth',1.5);
title('Kalman Filter 2 signals & sigma: 0.5')
legend('Sensor 1','Sensor 2','Combined Filtered Signal');

xlabel('Time (s)','fontsize',15);
ylabel('Voltage (V)','fontsize',15);



