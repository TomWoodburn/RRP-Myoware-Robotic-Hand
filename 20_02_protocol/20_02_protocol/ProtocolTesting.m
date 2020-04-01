%%
clear all;clc;close all
%% Light Tool Protocol

LT1 = importdata('20_02_Protocol_light_tool_1.csv');
LT2 = importdata('20_02_Protocol_light_tool_2.csv');

LTtime1=LT1.data(1417:end,1)';
LTvolt1=LT1.data(1417:end,2)';
LTVolt1N = normalize(LTvolt1,'scale');

LTtime2=LT2.data(:,1)';
LTvolt2=LT2.data(:,2)';
LTvolt2N = normalize(LTvolt2,'scale');

figure(1)
plot(LTtime1,LTVolt1N);grid on
hold on 

plot(LTtime2,LTvolt2N);grid on

hold off

%% Parallel Extension Protocol

PE1 = importdata('20_02_Protocol_parallel_extension_1.csv');
PE2 = importdata('20_02_Protocol_parallel_extension_2.csv');

PEtime1=PE1.data(:,1)';
PEvolt1=PE1.data(:,2)';
PEvolt1N = normalize(PEvolt1);

PEtime2=PE2.data(:,1)';
PEvolt2=PE2.data(:,2)';
PEvolt2N = normalize(PEvolt2);

figure(1)
plot(PEtime1,PEvolt1N);grid on
hold on 

plot(PEtime2,PEvolt2N);grid on

hold off

%% Medium Wrap Protocol (good signal)

MW1 = importdata('20_02_Protocol_Power_grasp_1_s.csv');
MW2 = importdata('20_02_Protocol_Power_grasp_2_s.csv');
MW3 = importdata('20_02_Protocol_Power_grasp_3.csv');

MWtime1=MW1.data(1072:2638,1)';
MWvolt1a=MW1.data(1072:2638,2)';
MWvolt1b=MW1.data(1072:2638,3)';
MWvolt1aN = normalize(MWvolt1a);
MWvolt1bN = normalize(MWvolt1b);

MWtime2=MW2.data(1072:2638,1)';
MWvolt2=MW2.data(1072:2638,2)';
MWvolt2N = normalize(MWvolt2);

MWtime3=MW3.data(1072:end,1)';
MWvolt3=MW3.data(1072:end,2)';
MWvolt3N = normalize(MWvolt3);


figure(1)


plot(MWtime1,MWvolt1a);grid on
hold on 

plot(MWtime1,MWvolt1b);grid on

%plot(MWtime3,MWvolt3);grid on

title('Selected Window: 2 Myoware Sensor Readings, "Medium Wrap" GRASP Taxonomy')
xlabel('Time from start of Protocol (Sec)')
ylabel('Recorded Voltage (V)')
legend('Signal 2: Centre of Palmarus Longus muscle','Signal 1: Middle of Extensor Digitorum muscle')

hold off

%% Medium Wrap Protocol 2 (dodgy signal)

MW1 = importdata('20_02_Protocol_Power_grasp_1_s.csv');
MW2 = importdata('20_02_Protocol_Power_grasp_2_s.csv');
MW3 = importdata('20_02_Protocol_Power_grasp_3.csv');

MWtime1=MW1.data(2638:end,1)';
MWvolt1a=MW1.data(2638:end,2)';
MWvolt1b=MW1.data(2638:end,3)';
MWvolt1aN = normalize(MWvolt1a);
MWvolt1bN = normalize(MWvolt1b);

MWtime2=MW2.data(2638:end,1)';
MWvolt2=MW2.data(2638:end,2)';
MWvolt2N = normalize(MWvolt2);

MWtime3=MW3.data(1072:end,1)';
MWvolt3=MW3.data(1072:end,2)';
MWvolt3N = normalize(MWvolt3);


figure(1)


plot(MWtime1,MWvolt1a);grid on
hold on 

plot(MWtime1,MWvolt1b);grid on

%plot(MWtime3,MWvolt3);grid on

title('Selected Window: 2 Myoware Sensor Readings, "Medium Wrap" GRASP Taxonomy')
xlabel('Time from start of Protocol (Sec)')
ylabel('Recorded Voltage (V)')
legend('Signal 2: Centre of Palmarus Longus muscle','Signal 1: Middle of Extensor Digitorum muscle')

hold off

%%
% 1. remove DC
min_s=min(s)
s=s+abs(min_s)
%%
% center signal, closest to dominant cycle to a sin cos signal
max_s=max(s)
s=s-max_s/2
hold on
plot(t,s);grid on;axis tight

%%
% 2. dominant base cycle w ith an fft:

S=fft(s)
absS=abs(S);
figure(2);
plot(absS);grid on;

max_absS=max(absS);
n_maxS=find(absS==max_absS)

n_maxS=n_maxS(1)
nT=floor(numel(t)/n_maxS)
%% Tom Edit using Squelch 

close all
%s=A(:,2)'
th1=0.3
figure(1);plot(t,s);grid on
[pks,locs]=findpeaks(s,'MinPeakHeight',th1)
hold on
plot(t(locs),pks,'bd')
hold off
%%
[pks,locs]=findpeaks(s,'MinPeakHeight',th1,'MinPeakDistance',floor(nT/2))
hold on
plot(t(locs),pks,'rd');

Freq = mean(diff(locs))

Var = var(diff(locs))^.5

dt=floor(.5*mean(diff(locs)))

sc2=zeros(numel(locs),2*dt+1)   % + - window span around each peak 
for k=1:1:numel(locs)
    if locs(k)<dt        % 1st cycle, with 1st peak closer than dt to beginning of signal.
       s0=s([1:locs(k)+dt])
       sc2(k,:)=[zeros(1,2*dt+1-numel(s0)) s0]
    end
      if locs(k)>numel(t)-dt  % last cycle with peak closer than dt to end of signal.
          sc2(k,:)=s([locs(k)-dt:end])
      end
      if locs(k)<numel(t)-dt && locs(k)>dt
          sc2(k,:)=s([locs(k)-dt:locs(k)+dt])
      end
end
% 4. and the eye diagram is:
figure(2);
for k=1:1:numel(locs)
    plot(sc2(k,:))
    hold on
end
grid on
% 005

%%
max_absS=max(absS);
n_maxS=find(absS==max_absS)

% 2nd value is just FFT mirror
n_maxS=n_maxS(1)   

max_amount_cycles=floor(numel(t)/2)

% if n_maxS = max_amount_cycles this would be the highest discernible
% frequency, with the FFT: it would be just 2 time samples per cycle.

% n_maxS: amount of cycles

hold on;plot(n_maxS,max_absS,'ro')

% nT: amount of samples per dominant cycle:
nT=floor(numel(t)/n_maxS)

% dominant cycle T
LTtime1=t([1:nT])
T=LTtime1(end)

% amount of lost samples ignoring .2

mod(nT,n_maxS)  % samples lost, no big deal

% instead of a for loop:
% sc=zeros(n_maxS,T)
% for ,,

LTvolt2=s([1:1:end-(numel(t)-n_maxS*nT)])


figure(3);
for k=1:1:3   % n_maxS
    plot(LTvolt2(k*[1:nT]))
    hold on;
end
% 003

% So the dominant cycle is not that constant.
% For the same basic cycle t(nT) we see 1,2 and 3 peaks
% This comes from a 2nd tone almost half way the dominant tone
% and really close to the dominant.
% fft approach is more reliable when there's a clear frequency peak
% but no high tones anywhere, particular near the dominant.

% the meaning of:
% repeating parts  peaks? valleys? up down transitions? down-up transitions?
% each containing a single action potential

% assuming the single action potential refers to a single peak



% 3. another way to get the basic cycle is with a squelch: you manually look at
% tjhe signal and set a threshold that includes all peaks above the
% threshold line. 

close all

s=LT1(:,2)'
th1=100
figure(1);plot(t,s);grid on

[pks,locs]=findpeaks(s,'MinPeakHeight',th1)
hold on
plot(t(locs),pks,'bd')


% avoid false peaks imposing min distance between found peaks, that may be
% for instance: nT

[pks,locs]=findpeaks(s,'MinPeakHeight',th1,'MinPeakDistance',floor(nT/2))
hold on
plot(t(locs),pks,'rd')
% 004

% So far what can we assert about this signal and the 'potential' events
% you have mentioned in the question?

% the events take place 
mean(diff(locs))

% fitter ~ standard deviation locs
var(diff(locs))^.5


dt=floor(.5*mean(diff(locs)))

sc2=zeros(numel(locs),2*dt+1)   % + - window span around each peak 
for k=1:1:numel(locs)
    if locs(k)<dt        % 1st cycle, with 1st peak closer than dt to beginning of signal.
       s0=s([1:locs(k)+dt])
       sc2(k,:)=[zeros(1,2*dt+1-numel(s0)) s0]
    end
    
    if locs(k)>numel(t)-dt  % last cycle with peak closer than dt to end of signal.
        sc2(k,:)=s([locs(k)-dt:end])
    end
    
    if locs(k)<numel(t)-dt && locs(k)>dt
        sc2(k,:)=s([locs(k)-dt:locs(k)+dt])
    end
        
end

% 4. and the eye diagram is:
figure(2);
for k=1:1:numel(locs)
    plot(sc2(k,:))
    hold on
end
grid on
% 005












% 6. comment, why the initial DC removal:
% The possible problem of using just 1 squelch level is that with noisy signals, most of real
% signals are, you set the squelch thinking it's going to work, you get a
% sample of N pulses ok, but next batch has a different noise ground level or a fluctuating
% noise ground, and there you are losing peaks.

% This is what happens when too much DC:
% sdc=A(:,2)'
% min_sdc=min(sdc)
% sdc=sdc+abs(min_sdc)
% Sdc=fft(sdc)
% absSdc=abs(Sdc);
% 
% figure(3);
% plot(absSdc);grid on;
% 
% max_absSdc=max(absSdc);
% n_maxSdc=find(absSdc==max_absSdc)
% 006

% n_maxSdc=1? no dominant cycle?? go for the 2nd? that may be confused for
% the 3rd? ..

