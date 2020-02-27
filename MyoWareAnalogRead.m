%% Acquire and analyze data from a temperature sensor

%% Connect to Arduino
% Use the arduino command to connect to an Arduino device.

a = arduino;

%% Take a single temperature measurement (V fixed)
% The datasheet for the TMP36 temperature sensor tells us that the voltage
% reading is directly proportional to temperature in Celsius with an 
% offset of 0.5V and a scale factor of 10 mV/°C (equivalent to 100 °C/V).
% Therefore the conversion can be represented as
%
% $T_C = (V-0.5)*100$
%
% We can read the output voltage, convert it to Celsius and convert the
% result to Farenheit as follows:
v = readVoltage(a,'A0');
%TempC = (v - 0.5)*100;
%TempF = 9/5*TempC + 32;
fprintf('VoltageReading from MyoWare Sensor:',v)

%% Record and plot 10 seconds of temperature data (V fixed)

ii = 0;
TempF = zeros(1e4,1);
t = zeros(1e4,1);

tic
while toc < 10
    ii = ii + 1;
    % Read current voltage value
    v(ii) = readVoltage(a,'A0');
    % Calculate temperature from voltage (based on data sheet)
    %TempC = (v - 0.5)*100;
    %TempF(ii) = 9/5*TempC + 32;
    % Get time since starting
    t(ii) = toc;
    %testing
end

% Post-process and plot the data. First remove any excess zeros on the
% logging variables.
v = v(1:ii);
t = t(1:ii);
% Plot temperature versus time
figure
plot(t,v,'-o')
xlabel('Elapsed time (sec)')
ylabel('Voltage (V)')
title('Myoware Voltage over time')
set(gca,'xlim',[t(1) t(ii)])

%% Compute acquisition rate

timeBetweenDataPoints = diff(t);
averageTimePerDataPoint = mean(timeBetweenDataPoints);
dataRateHz = 1/averageTimePerDataPoint;
fprintf('Acquired one data point per %.3f seconds (%.f Hz)\n',...
    averageTimePerDataPoint,dataRateHz)

%% Why is my data so choppy?

measurableIncrementV = 5/1023;
measurableIncrementC = measurableIncrementV*100;
measurableIncrementF = measurableIncrementC*9/5;
fprintf('The smallest measurable increment of this sensor by the Arduino is\n %-6.4f V\n %-6.2f°C\n %-6.2f°F\n',...
    measurableIncrementV,measurableIncrementC,measurableIncrementF);

%% Acquire and display live data (V fixed)

figure
hold on
h1 = animatedline();
h2 = animatedline('Color','r','LineWidth',0.5);
hold off
ax = gca;
ax.YGrid = 'on';
%ax.YLim = [65 85];

stop = false;
startTime = datetime('now');

while ~stop
    % Read current voltage value
    v1 = readVoltage(a,'A0');
    v2 = readVoltage(a, 'A1');
    % Calculate temperature from voltage (based on data sheet)
    %TempC = (v - 0.5)*100;
    %TempF = 9/5*TempC + 32;    
    % Get current time
    t =  datetime('now') - startTime;
    % Add points to animation
    addpoints(h1,datenum(t),v1)
    addpoints(h2,datenum(t),v2)
    % Update axes
    ax.XLim = datenum([t-seconds(15) t]);
    datetick('x','keeplimits')
    drawnow
    % Check stop condition
    stop = readDigitalPin(a,'D12');
end

%% Plot the recorded data 

[timeLogs1,tempLogs1] = getpoints(h1);
[timeLogs2,tempLogs2] = getpoints(h2);
timeSecs1 = (timeLogs1-timeLogs1(1))*24*3600;
timeSecs2 = (timeLogs2-timeLogs2(1))*24*3600;
figure
hold on
plot(timeSecs1,tempLogs1)
plot(timeSecs2,tempLogs2)
hold off
xlabel('Elapsed time (sec)')
ylabel('Recorded Voltage (V)')

%% Smooth out readings with moving average filter

smoothedTemp = smooth(tempLogs1,25);
tempMax = smoothedTemp + 2*9/5;
tempMin = smoothedTemp - 2*9/5;

figure
plot(timeSecs,tempLogs1, timeSecs,tempMax,'r--',timeSecs,tempMin,'r--')
xlabel('Elapsed time (sec)')
ylabel('Recorded Voltage (V)')
hold on 

%%
% Plot the original and the smoothed temperature signal, and illustrate the
% uncertainty.

plot(timeSecs,smoothedTemp,'r')

%% Save results to a file

T = table(timeSecs1',tempLogs2',tempLogs1','VariableNames',{'Time_sec','Voltage1', 'Voltage2'});
filename = 'Protocol_20_02_parallel_extension_2_sensors_2_switched.xlsx';
% Write table to file 
writetable(T,filename)
% Print confirmation to command line
fprintf('Results table with %g temperature measurements saved to file %s\n',...
    length(timeSecs1),filename)
