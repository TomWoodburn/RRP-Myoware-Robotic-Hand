clear all; clc; close all;

% Discretization of continuous system
dt = 0.01;
data = csvread('20_02_protocol/20_02_Protocol_light_tool_1.csv',1,0);
% Generate time array
data = data(1400: size(data),:);
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
Sigma = 0.8;
R=[(Sigma^2) 0;
    0 (0.5^2)];

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

splits = 60;

itp = findchangepts(x(1,:),'MaxNumChanges',splits,'Statistic','mean');


figure();set(gcf,'color','white');
hold on;
itptrimpre = [1];
for i= 1:splits-1
    if itp(i+1)-itp(i) >60
        itptrimpre = [itptrimpre itp(i+1)];
    end
    
end
itptrim=[];
for i = 1:size(itptrimpre,2)-1
    if mean(x(1,itptrimpre(i):itptrimpre(i+1))) > 0.5
        itptrim = [itptrim [itptrimpre(i);1]];
    else
        itptrim = [itptrim [itptrimpre(i);0]];
    end
end
itp = itptrim;
plot(t,SignalNoisy(1,:),'r');
plot(t,SignalNoisy(2,:));
itpSize = size(itp,2);
for i=1:itpSize
    if itp(2,i)==1
        xline(t(itp(1,i)));
    else
        xline(t(itp(1,i)),'r');
    end
end
plot(t,x(1,:))
hold off

%%
for i =1:itpSize-1
    figure
    hold on
    %lineGraph
    %plot(SignalNoisy(1,itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)))
    %ScatterGraph
    scatter(SignalNoisy(1,itp(1,i):itp(1,i+1)),SignalNoisy(2,itp(1,i):itp(1,i+1)))
    hold off
    %The splits
    figure
    hold on
    plot(t(itp(1,i):itp(1,i+1)),SignalNoisy(1,itp(1,i):itp(1,i+1)),'r');
    plot(t(itp(1,i):itp(1,i+1)),SignalNoisy(2,itp(1,i):itp(1,i+1)));
    X = [SignalNoisy(1,itp(1,i):itp(1,i+1));SignalNoisy(2,itp(1,i):itp(1,i+1))]';
    %PolyFit
    j=6
    
    [p, S] = polyfit(t(itp(1,i):itp(1,i+1)),SignalNoisy(2,itp(1,i):itp(1,i+1)),j);
    f = polyval(p,t(itp(1,i):itp(1,i+1)));
    plot(t(itp(1,i):itp(1,i+1)),SignalNoisy(2,itp(1,i):itp(1,i+1)),'*',t(itp(1,i):itp(1,i+1)),f,'-')
    title(['Order is: ',num2str(j),' norm is: ',num2str(S.normr)])

    hold off
    
end
