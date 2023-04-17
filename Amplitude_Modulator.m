clf
%Graph Parameters
periods = 3 %Will visualize 3 periods of the message frequency
resolution = 0.0000001 %This is the resolution of the wave. Smaller number = better resolution. Recommended at least 0.0001

%Inputs to the equations
carrierdBm = 0;     %dBm
DCOffset = 10.167; %Volts
R1 = 10000
R2 = 1200
messageFrequency = 1000000;       %The frequency for the message is 1MHz
carrierFrequency = 10000000000; %The frequency for the carrier wave is 10GHz
impedance = 50;     %Ohms
transferFuncDataPoints = [-1.34 -1.22 -1.1;
                          20    15   9]; %A 2d array with datapoints about the transfer function.
             % Voltage  Gain(dB)  Used to calculate the transfer function
             %To make the transfer function more accurate, take more
             %measurements of the graph "Gain vs Control Voltage" on
             %datasheet


%Formulas%
T=1/messageFrequency; %Obtains the period of the message frequency
t = 0:resolution*T:periods*T; %Creates the vector that represents the time
carrierWave = ((sqrt(impedance/1000)*10^(carrierdBm/20)*cos(2*pi*carrierFrequency*t))*2)/sqrt(2); %The carrier wave as Vpeak
messageWave = .5*cos(2*pi*messageFrequency*t) + .5*cos(4*pi*messageFrequency*t); %Tonal Message


%The carrier Frequency vs time
plot(t, carrierWave);
title("Carrier Wave (10GHz @ 5 dBm)")
xlabel("Time (Seconds)");
ylabel("Voltage (V)");
ax = gcf;
exportgraphics(ax, "CarriervTime.png");

%The message vs time
plot(t, messageWave);
title("Message Wave (Tonal Frequency)")
xlabel("Time (Seconds)");
ylabel("Voltage (V)");
ax = gcf;
exportgraphics(ax, "MessagevTime.png");

%Peform Linear Regression using the noted date points to create the
%transfer function
regressionX = transferFuncDataPoints(1,:);
regressionY = transferFuncDataPoints(2,:);
mdl = fitlm(regressionX, regressionY);

title("Linear Region: Gain vs Control Voltage (Transfer Function) HMC694LP4")
xlabel("Control Voltage [V]");
ylabel("Gain [dB]");
ax = gcf;
exportgraphics(ax, "LinearRegion.png");

%Extracted coefficients in the form y=mx+b
m = mdl.Coefficients{2, 1}; %Has units of dB/V
b = mdl.Coefficients{1, 1}; %Has units of dB
fprintf("The regression is calcualted to be y = %fx + %f\n", m, b);

%Plot the fit of the linear regression
scatter(regressionX, regressionY);
hold on
fitX = [regressionX(1, 1):0.01:regressionX(1, end)];
%dataPoints = -1.34:0.01:-1.1
fitY = m*fitX + b
plot(fitX, fitY)
title("Fitting the Linear Region of the Variable Gain Amplifier")
xlabel("Control Voltage (V)");
ylabel("Gain (dB)");
ax = gcf;
exportgraphics(ax, "LinearRegression.png");

hold off
%Now, let's attenuate the Message Frequency and DC offset
%This is representative of the inverting
%summing amplifier, which will sum
%messageWave = (-.22/1.5625) * messageWave - 1.2108
messageWave = -(R2/R1)*((messageWave)+(DCOffset))
plot(t, messageWave);
title("Message after Inverting Summing Amplifier")
xlabel("Time (Seconds)");
ylabel("Voltage (V)");
ax = gcf;
exportgraphics(ax, "AfterInverting.png");

%The linear regression transfer function represents the Analog Variable
%Gain Amplifier. The regression will match the control voltage to its
%corresponding gain (dB). Thus, we will pass the equation into the transfer
%function.
transferFunctionGain = m * messageWave + b;
plot(t, transferFunctionGain);
title("Transfer Function in the Linear Region")
xlabel("Control Voltage (V)");
ylabel("Gain (dB)");
ax = gcf;
exportgraphics(ax, "Message Gain.png");
%Once we have that, let's add it to the corresponding wave from the local
%oscillator. I'm going to transfer everything into dBms so that we can add
%it up. 

AMWaveV = carrierWave .* 10.^(transferFunctionGain/20);
plot(t, AMWaveV);
title("AM Wave")
xlabel("Time (Seconds)");
ylabel("Voltage (V)");
ax = gcf;
exportgraphics(ax, "AMWave.png");

%This is still centered at 15dB, not 15dBm
%Therefore, we attenuate it with 30dBm, giving us this equation

FinalWave = AMWaveV ./10.^(30/20);
plot(t,FinalWave);
title("Final Wave")
xlabel("Time (Seconds)");
ylabel("Voltage (V)");
ax = gcf;
exportgraphics(ax, "FinalWave.png");