%% Plotting, Averaging & Sampling the Protocol Signals. 27/02/2020
%% Plotting the CSVs

%Array=csvread('20_02_Protocol_light_tool_1.csv');

LT1 = importdata('20_02_protocol/20_02_Protocol_light_tool_1.csv');
col1 = LT1.data(1417:end, 1)';
col2 = LT1.data(1417:end, 2)';
col3 = LT1.data(1417:end, 3)';

%Normcol2 = normalize(col2, 'norm', 1);

%Normcol2 = col2 ./ max(col2,'includenan','all')

Normalize(col2,'range');

%% Plotting Gradient of the curve
FY = gradient(col2);

figure()
[TF, S1] = ischange(FY);
plot(FY,'*')
hold on
plot(col2)
hold on
stairs(S1)
legend('Data','Segment Mean','Location','NW')
hold off


%method = mean, variance, linear

figure()
%[TF, S2, S3] = ischange(col2,'mean','theshold',200);
[TF,S2,S3] = ischange(col2,'linear','MaxNumChanges',10);
plot(TF)
hold on
plot(col2)
hold on
stairs(S2)
legend('Data','Segment Mean','Location','NW')
hold off



%%
%plot(col1, col2)
figure(2)
tiledlayout(2,1)
nexttile
plot(col1,col2,'r-','linewidth',3)
hold on;
plot(col1,col3,'b-','linewidth',3)
xlabel('Sample','fontsize',16)
ylabel('Voltage from Myoware Sensor (V)','fontsize',16)
set(gca,'fontsize',16)
hold off

nexttile
[SignChange] = ischange(FY,'linear','Threshold',0.01);
%segline = S1.*(1:3197) + S2;
plot(SignChange)              
legend('Data')

%% Change in Linear Regime

[SignChange] = ischange(FY,'linear','Threshold',0.01);
%segline = S1.*(1:3197) + S2;
plot(SignChange)              
legend('Data')

%% findchangepts

figure()
%title('Detecting Changes RMS')
findchangepts(col2,'MaxNumChanges',34,'Statistic','rms')

figure()
%title('Detecting Changes Mean')
findchangepts(col2,'MaxNumChanges',34,'Statistic','mean')


figure()
%title('Detecting Changes std')
findchangepts(col2,'MaxNumChanges',34,'Statistic','std')


figure()
%title('Detecting Changes linear')
findchangepts(col2,'MaxNumChanges',34,'Statistic','linear')

%}

%% IsChange - Finding Abrupt Changes in Data 

