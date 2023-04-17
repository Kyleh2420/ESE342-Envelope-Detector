%Global Variables
inputFile = "Gonna.mp3";
outputFile = "output.wav";
modulationIndex = 0.3;
phaseAngle = 0;
carrierFrequency = 560*10^3; %560kHz
R = 75000;
C = 1*10^-9;
interpolatedSample = 15;


useSample = false;
sampleFrequency = 440;
sampleR = 700;
sampleC = 1*10^-6;


%Read the audio file. Then, transpose and convert to single channel.
%y: Matrix containing the audio data. See: https://www.mathworks.com/help/matlab/ref/audioread.html#btiabil-1-y
%For MP3s: Audio data ranges between -1 and 1. 
%Fs: Sample Rate
[y, Fs] = audioread(inputFile);
y = y.';
y = mean(y);
t = 0:1/Fs:(length(y)/Fs)-(1/Fs);  

if useSample 
    y = sin(2*pi*sampleFrequency*t);
    z=y;
    R = sampleR;
    C = sampleC;
    sound(y, Fs);
end


plot(t, y);
title("Audio file Intensity vs Time")
xlabel("Time (seconds)");
ylabel("Internsity");
if useSample
    axis([0 1/sampleFrequency -1 1]);
end
ax = gcf;
exportgraphics(ax,"Audio V Time.png");
%Creates the vector that represents the time at each sample, from 0 to the
%end

%Interpolation of Data
tInterpolated = 0:1/(interpolatedSample*carrierFrequency):(length(y)/Fs)-(1/Fs);
yInterpolated = interp1(t, y, tInterpolated);


plot(tInterpolated, yInterpolated);
title("Audio file Interpolation vs Time")
xlabel("Time (seconds)");
ylabel("Internsity");
if useSample
    axis([0 1/sampleFrequency -1 1])
end
ax = gcf;
exportgraphics(ax,"Audio V Time Interpolated.png");

%Begin Attempted FFT section

%{
%Sampling Frequency is 1/(interpolatedSample*carrierFrequency)
NFFT = length(yInterpolated);
%Length of signal is the size of the arrays
audioFFT = fft(yInterpolated, NFFT);
F = ((0:1/NFFT:1-1/NFFT)*(1/(interpolatedSample*carrierFrequency))).';

magnitudeY = abs(audioFFT);        % Magnitude of the FFT
phaseY = unwrap(angle(audioFFT));  % Phase of the FFT

%helperFrequencyAnalysisPlot1(F,magnitudeY,phaseY,length(yInterpolated));
dB_mag=mag2db(magnitudeY);
subplot(2,1,1);plot(F(1:NFFT/2),dB_mag(1:NFFT/2));title('Magnitude response of signal');
ylabel('Magnitude(dB)');
subplot(2,1,2);plot(F(1:NFFT/2),phaseY(1:NFFT/2));title('Phase response of signal');
xlabel('Frequency in kHz')
ylabel('radians');

cla()

%}
%End Attempted FFT section

vc = 5.0*(1+(modulationIndex.*yInterpolated)).*cos(2*pi*carrierFrequency.*tInterpolated + phaseAngle)*(10^-6);

plot(tInterpolated, vc);
title("Modulated Signal Voltage vs Time")
xlabel("Time (seconds)");
ylabel("Voltage (V)");
ax = gcf;
if useSample
    axis([0 1/sampleFrequency -inf inf])
end
exportgraphics(ax,"Modulated Signal Before Normalization.png");

%I don't think we can do this. If we multiply the x(t) by anything, it'll
%screw with mu, the modulation index. However, what we do need to do is
%normalize x(t) to go from 1 to -1
%yInterpolated = yInterpolated*100000;

yInterpolated = yInterpolated * (1/max(yInterpolated));

%Creating the modulated signal, units are V

vc = 5.0*(1+(modulationIndex.*yInterpolated)).*cos(2*pi*carrierFrequency.*tInterpolated + phaseAngle)*(10^-6);

plot(tInterpolated, vc);
title("Modulated Signal Voltage vs Time")
xlabel("Time (seconds)");
ylabel("Voltage (V)");
ax = gcf;
if useSample
    axis([0 1/sampleFrequency -inf inf])
end
exportgraphics(ax,"Modulated Signal After Normalization.png");

%Now that the signal is transmitted over the waves, lets amplify at the
%reciever module by a factor of 10^6

vc = vc*100000;

%Add the DC Bias, 0.7V
bias = vc + 0.7;
plot (tInterpolated, bias);
title("After 1V Bias voltage has been added")
xlabel("Time (seconds)");
ylabel("Voltage (V)");
ax = gcf;
if useSample
    axis([0 1/sampleFrequency -1 1])
end
exportgraphics(ax,"Bias.png");


diode = zeros(1,length(bias));
%Rectify the signal using an ideal practical diode with a forward drop of 1V. 
diode = bias - 0.7;
for c = 1:length(bias)
    if diode(1,c) <= 0
        diode(1,c) = 0;
    end
end

%t = linspace(0,size(diode,2),size(diode,2));
%subplot(2,2,3)
plot(tInterpolated, diode);
title("Diode vs Time")
xlabel("Time (seconds)");
ylabel("Voltage (V)");
ax = gcf;
if useSample
    axis([0 1/sampleFrequency -1 1])
    %axis([0 0.0001 0 0.15])
end
exportgraphics(ax, "Diode Signal.png");
%axis([0 inf 0  0.000006])


%Then, we pass it through a simple peak detector consisting of a resistor
%and capacitor in parallel. The discharge of the capacitor is given to be
%V(t) = V0*e^(-t/(RC)), and the upcharge will follow the peaks. Thus, a
%piecewise function is formed.
%Charging the cap (Or when V(n-1) < V(n))
%Discharging the cap (Or when V(n-1) >= V(n))
peakDetector = zeros(1,length(diode));
peakDetector(1, 1) = diode (1, 1);
V0 = 0;
recentlyCharged = 0;
timeSinceCharge = 0;
activeCapacitor = 0;
for c = 2:length(diode)
    if diode(1, c-1) < diode(1, c) && activeCapacitor == 0
        peakDetector(1, c) = diode(1, c);
        V0 = diode(1, c);
        recentlyCharged = c;
    else
        %To find t, we need to figure out the number of samples since the
        %last charge of the capacitor, then divide by the samples per
        %second to get the time.
        %When the active capacitor discharges, we have to wait until it hits the function again
        timeSinceCharge = c - recentlyCharged; %Measured in Samples
        timeSinceCharge = timeSinceCharge/Fs; %Divided by samples per second to get seconds
        peakDetector(1, c) = V0 * exp(-(timeSinceCharge)/(R*C));
        activeCapacitor = 1;
        if peakDetector(1, c) <= diode(1, c)
            activeCapacitor = 0;
        end
    end
end
%Change the time vector

%t = linspace(0,size(peakDetector,2),size(peakDetector,2));
%subplot(2,2,4)
plot(tInterpolated, peakDetector);
title("Demodulated Signal Voltage vs Time")
xlabel("Time (seconds)");
ylabel("Voltage (V)");
ax = gcf;
if useSample
    axis([0 1/sampleFrequency -1 1])
    %axis([0.0016 0.002 0 0.15])
end
exportgraphics(ax, "Demodulated Signal.png");



x = round(carrierFrequency * interpolatedSample/Fs);
convolutionTable = zeros(1, x);
for c = 1:length(convolutionTable)
    convolutionTable(1,c) = 1/x;
end

counter = 0;
j = 1;
avg = 0;
output = zeros(1, 1);
for i = 1:length(peakDetector)
    avg = avg + peakDetector(1, i);
    counter = counter + 1;
    if (counter == x)
        output(1, j) = avg/x;
        j = j+1;
        counter = 0;
        avg = 0;
    end
end
%I tried to do it the smart way
%Add 0s until we can divide by x to get a whole number
%closeToEnd = mod(peakDetectSize, x);
%for i = size(peakDetector,2)-closeToEnd:size(peakDetector,2)+1
%    peakDetector(1, i) = 0;
    
%end

%audioAverage = reshape(peakDetector, round(size(peakDetector,2)/x), x);
%pspectrum(peakDetector,10000,"spectrogram")


sound(output, Fs);
audiowrite(outputFile, output, Fs);