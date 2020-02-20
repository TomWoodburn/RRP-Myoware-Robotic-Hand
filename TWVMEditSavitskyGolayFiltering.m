%Savitsky-Golay filtering
%%
figure(1);
x = importdata('Protocol1Edit.csv');
th = x.data(:,2);
time = x.data(:,1);
plot(time,th,'k-','linewidth',3);
xlabel('Time (Sec)','fontsize',16)
ylabel('Voltage from Myoware [V]','fontsize',16)
set(gca,'fontsize',16)

%%
thf = sgolayfilt(th,7,119); % Savitsky-Golay filter
figure(4);
plot(time,th,'k-','linewidth',0.5);
hold on;
plot(time,thf,'g-','linewidth',0.5);
legend('Raw data','Filtered data');
xlabel('Time (Sec)','fontsize',16)
ylabel('Voltage from Myoware [V]','fontsize',16)
set(gca,'fontsize',16)
%hold off;