%%Simple smoothing
figure(1);
x = importdata('Protocol1Edit.csv');
th = x.data(:,2);
time = x.data(:,1);
plot(time,th,'k-','linewidth',0.5);
xlabel('Time (Sec)','fontsize',16)
ylabel('Voltage from Myoware [V]','fontsize',16)
set(gca,'fontsize',16)
%Then move a window while taking average in that window
win = [-15:15]
sth = [th(1:30)];
for i = 16:length(th)-16
  sth(i) = mean(th(i+win));  
end
figure(2);
plot(time,th,'b-','linewidth',0.5);
hold on;
plot(sth,'r-','linewidth',0.5);
legend('Raw data','Smoothed data');
xlabel('Time (Sec)','fontsize',16)
ylabel('Voltage from Myoware [V]','fontsize',16)
set(gca,'fontsize',16)

% In Matlab, smoothing can be done using a single command
% Z = smooth(Y,SPAN)
sth = smooth(th,30);
figure(3);
plot(th,'b-','linewidth',0.5);
hold on;
plot(sth,'r-','linewidth',0.5);
legend('Raw data','Smoothed data');
xlabel('Time (Sec)','fontsize',16)
ylabel('Voltage from Myoware [V]','fontsize',16)
set(gca,'fontsize',16)