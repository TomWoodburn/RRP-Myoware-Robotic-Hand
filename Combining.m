clear all; clc; close all;

% Discretization of continuous system
dt = 0.01;
dataLT = csvread('20_02_protocol/20_02_Protocol_light_tool_1.csv',1,0);
% Generate time array
dataLT = dataLT(1400: size(dataLT),:);
tlt = dataLT(:,1)';

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
SignalNoisy = dataLT(:,2:3)';

% LQE function
xl = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);


X = [sgolayfilt(SignalNoisy(1,:),7,111); sgolayfilt(SignalNoisy(2,:),7,111)]';
GMModel = fitgmdist(X,2);
figure('Name','Light tool')
y = [zeros(1000,1);ones(1000,1)];
h = gscatter(X(:,1),X(:,2));
hold on
gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
axis([0 inf 0 5]);
hold off
%% Plots


figure();set(gcf,'color','white');
hold on;

plot(tlt,SignalNoisy(1,:),'r');
plot(tlt,SignalNoisy(2,:));
plot(tlt,xl(1,:),'b','linewidth',1.5);



splitsL = 60;
itpL = findchangepts(xl(1,:),'MaxNumChanges',splitsL,'Statistic','mean');
figure();set(gcf,'color','white');
hold on;
itptrimpre = [1];
for i= 1:splitsL-1
    if itpL(i+1)-itpL(i) >60
        itptrimpre = [itptrimpre itpL(i+1)];
    end
    
end
itptrim=[]
for i = 1:size(itptrimpre,2)-1
    if mean(xl(1,itptrimpre(i):itptrimpre(i+1))) > 0.5
        itptrim = [itptrim [itptrimpre(i);1]];
    else
        itptrim = [itptrim [itptrimpre(i);0]];
    end
end
itpL = itptrim;
plot(tlt,SignalNoisy(1,:),'r');
plot(tlt,SignalNoisy(2,:));
itpSize = size(itpL,2);
for i=1:itpSize
    if itpL(2,i)==1
        xline(tlt(itpL(1,i)));
    else
        xline(tlt(itpL(1,i)),'r');
    end
end
plot(tlt,xl(1,:))

hold off
%%

for i =1:itpSize-1
    figure
    hold on
    %lineGraph
    %plot(SignalNoisy(1,itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)))
    %ScatterGraph
    scatter(SignalNoisy(1,itpL(1,i):itpL(1,i+1)),SignalNoisy(2,itpL(1,i):itpL(1,i+1)))
    hold off
    %The splits
    figure
    hold on
    plot(tlt(itpL(1,i):itpL(1,i+1)),SignalNoisy(1,itpL(1,i):itpL(1,i+1)),'r');
    plot(tlt(itpL(1,i):itpL(1,i+1)),SignalNoisy(2,itpL(1,i):itpL(1,i+1)));
    X = [sgolayfilt(SignalNoisy(1,itpL(1,i):itpL(1,i+1)),5,7);sgolayfilt(SignalNoisy(2,itpL(1,i):itpL(1,i+1)),5,7)]';
    %PolyFit
    %[p, S, mu] = polyfit(t(itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)),10);
    %plot(t(itp(i):itp(i+1)), polyval(p, t(itp(i):itp(i+1)))); 
    hold off
    GMModel = fitgmdist(X,1);
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
end
%%
dataPG = csvread('20_02_protocol/20_02_Protocol_Power_grasp_3.csv',1,0);
dataPG = dataPG(2300:5804,:);
SignalNoisy = dataPG(:,2:3)';
tpg = dataPG(:,1)';
xp = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);
X = [sgolayfilt(SignalNoisy(1,:),7,111); sgolayfilt(SignalNoisy(2,:),7,111)]';
GMModel = fitgmdist(X,2);
figure('Name','Power grasp')
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
figure();set(gcf,'color','white');
hold on
plot(tpg,SignalNoisy(1,:),'r');
plot(tpg,SignalNoisy(2,:));
plot(tpg,xp(1,:),'b','linewidth',1.5);
splitsPG = 60;
itpPG = findchangepts(xp(1,:),'MaxNumChanges',splitsPG,'Statistic','mean');
figure();set(gcf,'color','white');
hold on;

itptrimpre = [1];
for i= 1:splitsPG-1
    if itpPG(i+1)-itpPG(i) >60
        itptrimpre = [itptrimpre itpPG(i+1)];
    end
    
end
itptrim=[]
for i = 1:size(itptrimpre,2)-1
    if mean(xp(1,itptrimpre(i):itptrimpre(i+1))) > 0.5
        itptrim = [itptrim [itptrimpre(i);1]];
    else
        itptrim = [itptrim [itptrimpre(i);0]];
    end
end


itpPG = itptrim;
plot(tpg,SignalNoisy(1,:),'r');
plot(tpg,SignalNoisy(2,:));
itpSize = size(itpPG,2);
for i=1:itpSize
    if itpPG(2,i)==1
        xline(tpg(itpPG(1,i)));
    else
        xline(tpg(itpPG(1,i)),'r');
    end
end
plot(tpg,xp(1,:))
hold off
%%

for i =1:itpSize-1
    figure
    hold on
    %lineGraph
    %plot(SignalNoisy(1,itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)))
    %ScatterGraph
    scatter(SignalNoisy(1,itpPG(1,i):itpPG(1,i+1)),SignalNoisy(2,itpPG(1,i):itpPG(1,i+1)))
    hold off
    %The splits
    figure
    hold on
    plot(tpg(itpPG(1,i):itpPG(1,i+1)),SignalNoisy(1,itpPG(1,i):itpPG(1,i+1)),'r');
    plot(tpg(itpPG(1,i):itpPG(1,i+1)),SignalNoisy(2,itpPG(1,i):itpPG(1,i+1)));
    X = [sgolayfilt(SignalNoisy(1,itpPG(1,i):itpPG(1,i+1)),5,7);sgolayfilt(SignalNoisy(2,itpPG(1,i):itpPG(1,i+1)),5,7)]';
    %PolyFit
    %[p, S, mu] = polyfit(t(itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)),10);
    %plot(t(itp(i):itp(i+1)), polyval(p, t(itp(i):itp(i+1)))); 
    hold off
    GMModel = fitgmdist(X,1);
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
end
%%
dataS = csvread('20_02_protocol/20_02_Protocol_sphere_finger_1.csv',1,0);
% Generate time array
dataS = dataS(800:4700,:);
ts = dataS(:,1)';
SignalNoisy = dataS(:,2:3)';
xs = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);
X = [sgolayfilt(SignalNoisy(1,:),7,111); sgolayfilt(SignalNoisy(2,:),7,111)]';
GMModel = fitgmdist(X,2);
figure('Name','Sphere finger')
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
figure();set(gcf,'color','white');
hold on
plot(ts,SignalNoisy(1,:),'r');
plot(ts,SignalNoisy(2,:));
plot(ts,xs(1,:),'b','linewidth',1.5);
splitsS = 60;
itpS = findchangepts(xs(1,:),'MaxNumChanges',splitsS,'Statistic','mean');
figure();set(gcf,'color','white');
hold on;
itptrimpre = [1];
for i= 1:splitsS-1
    if itpS(i+1)-itpS(i) >60
        itptrimpre = [itptrimpre itpS(i+1)];
    end
    
end
itptrim=[]
for i = 1:size(itptrimpre,2)-1
    if mean(xs(1,itptrimpre(i):itptrimpre(i+1))) > 0.5
        itptrim = [itptrim [itptrimpre(i);1]];
    else
        itptrim = [itptrim [itptrimpre(i);0]];
    end
end


itpS = itptrim;
plot(ts,SignalNoisy(1,:),'r');
plot(ts,SignalNoisy(2,:));
itpSize = size(itpS,2);
for i=1:itpSize
    if itpS(2,i)==1
        xline(ts(itpS(1,i)));
    else
        xline(ts(itpS(1,i)),'r');
    end
end
plot(ts,xs(1,:))
hold off
%%

for i =1:itpSize-1
    figure
    hold on
    %lineGraph
    %plot(SignalNoisy(1,itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)))
    %ScatterGraph
    scatter(SignalNoisy(1,itpS(1,i):itpS(1,i+1)),SignalNoisy(2,itpS(1,i):itpS(1,i+1)))
    hold off
    %The splits
    figure
    hold on
    plot(ts(itpS(1,i):itpS(1,i+1)),SignalNoisy(1,itpS(1,i):itpS(1,i+1)),'r');
    plot(ts(itpS(1,i):itpS(1,i+1)),SignalNoisy(2,itpS(1,i):itpS(1,i+1)));
    X = [sgolayfilt(SignalNoisy(1,itpS(1,i):itpS(1,i+1)),5,9);sgolayfilt(SignalNoisy(2,itpS(1,i):itpS(1,i+1)),5,9)]';
    %PolyFit
    %[p, S, mu] = polyfit(t(itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)),10);
    %plot(t(itp(i):itp(i+1)), polyval(p, t(itp(i):itp(i+1)))); 
    hold off
    GMModel = fitgmdist(X,1);
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
end
%%
dataPE = csvread('20_02_protocol/20_02_Protocol_parallel_extension_1.csv',1,0);
% Generate time array
dataPE = dataPE(1050:4400,:);
tpe = dataPE(:,1)';
SignalNoisy = dataPE(:,2:3)';
xpe = KalmanFilter(A,C,Q,R,xInit,PInit,SignalNoisy);
X = [sgolayfilt(SignalNoisy(1,:),7,111); sgolayfilt(SignalNoisy(2,:),7,111)]';
GMModel = fitgmdist(X,2);
figure('Name','Parallel extension')
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
figure();set(gcf,'color','white');
hold on
plot(tpe,SignalNoisy(1,:),'r');
plot(tpe,SignalNoisy(2,:));
plot(tpe,xpe(1,:),'b','linewidth',1.5);
splitsPE = 90;
itpPE = findchangepts(xpe(1,:),'MaxNumChanges',splitsPE,'Statistic','mean');
figure();set(gcf,'color','white');
hold on;

itptrimpre = [1];
for i= 1:splitsPE-1
    if itpPE(i+1)-itpPE(i) >50
        itptrimpre = [itptrimpre itpPE(i+1)];
    end
    
end

itptrim=[]
for i = 1:size(itptrimpre,2)-1
    if mean(xpe(1,itptrimpre(i):itptrimpre(i+1))) > 0.5
        itptrim = [itptrim [itptrimpre(i);1]];
    else
        itptrim = [itptrim [itptrimpre(i);0]];
    end
end

itpPE = itptrim;
plot(tpe,SignalNoisy(1,:),'r');
plot(tpe,SignalNoisy(2,:));
itpSize = size(itpPE,2);
for i=1:itpSize
    if itpPE(2,i)==1
        xline(tpe(itpPE(1,i)));
    else
        xline(tpe(itpPE(1,i)),'r');
    end
end
plot(tpe,xpe(1,:))
hold off

%%

for i =1:itpSize-1
    figure
    hold on
    %lineGraph
    %plot(SignalNoisy(1,itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)))
    %ScatterGraph
    scatter(SignalNoisy(1,itpPE(1,i):itpPE(1,i+1)),SignalNoisy(2,itpPE(1,i):itpPE(1,i+1)))
    hold off
    %The splits
    figure
    hold on
    plot(tpe(itpPE(1,i):itpPE(1,i+1)),SignalNoisy(1,itpPE(1,i):itpPE(1,i+1)),'r');
    plot(tpe(itpPE(1,i):itpPE(1,i+1)),SignalNoisy(2,itpPE(1,i):itpPE(1,i+1)));
    X = [sgolayfilt(SignalNoisy(1,itpPE(1,i):itpPE(1,i+1)),5,7);sgolayfilt(SignalNoisy(2,itpPE(1,i):itpPE(1,i+1)),5,7)]';
    %PolyFit
    %[p, S, mu] = polyfit(t(itp(i):itp(i+1)),SignalNoisy(2,itp(i):itp(i+1)),10);
    %plot(t(itp(i):itp(i+1)), polyval(p, t(itp(i):itp(i+1)))); 
    hold off
    GMModel = fitgmdist(X,1);
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
end
%% Combining
data = [dataLT; dataPG; dataPE; dataS];
itpPE = itpPE + [size(dataLT,1);0];
itpPG = itpPG + [size(dataLT,1);0]+[size(dataPE,1);0];
itpS = itpS + [size(dataLT,1);0]+[size(dataPE,1);0] + [size(dataPG,1);0];

tpe = tpe+tpe(size(dataLT,1));
tpg = tpg+tlt(size(dataLT,1))+tpe(size(dataPE,1));
ts = ts+tlt(size(dataLT,1))+tpe(size(dataPE,1))+tpg(size(dataPG,1));

SignalNoisy = data(:,2:3)';
t = [tlt tpe tpg ts];
figure()
hold on
plot(t,SignalNoisy(1,:),'r');
plot(t,SignalNoisy(2,:));
itp = [itpL itpPE itpPG itpS];
%%
figure
hold on
trueData = [];
truet = []
for i = 1:size(itp,2)
    if itp(2,i) == 1
        trueData = [trueData; data(itp(1,i):itp(1,i+1),:)];
        truet = [truet,  t(itp(1,i):itp(1,i+1))];
    end
end 
plot(truet,trueData(:,3))
plot(truet,trueData(:,2))
hold off
%%
SignalNoisy = sgolayfilt(data(:,2:3),5,149);
X = SignalNoisy;
GMModel = fitgmdist(X,4);
figure()
y = [zeros(1000,1);ones(1000,1)];
h = gscatter(X(:,1),X(:,2));
hold on
gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
axis([0 inf 0 5])
hold off