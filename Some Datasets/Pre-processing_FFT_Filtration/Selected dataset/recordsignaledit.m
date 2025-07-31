% This code config the Fireface device for sending and reciving data using DSP library (digital signal processing).
% Realizing essential parameter for data transfer https://www.mathworks.com/help/dsp/ug/sample-and-frame-based-concepts.html

% close all
clear all
clc
% config parameters
DeviceName = 'ASIO Fireface USB';
dtyp='32-bit float';
fs = 96000; % Available sample rates: 32000, 44100, 48000, 64000, 88200, 96000, 128000, 176400, 192000
input_channel = 5; % Available setting is 5 6 7 8
BufferSize  = 4096; % Max is 4096

frameSize = 1024;  %2048 make sure that this parameter is high enough for sending data in diff. freq

% Configuration for Recording signal from the ME sensor
recorder = dsp.AudioRecorder('SamplesPerFrame', frameSize, 'DeviceDataType','32-bit float', 'SampleRate',fs);
recorder.DeviceName =  'ASIO Fireface USB';
recorder.NumChannels = input_channel; 
                         
% Configuration for Sending the excitation signal to the coil
setpref('dsp','portaudioHostApi',3)  % Output channel more info https://www.mathworks.com/help/dsp/ref/toaudiodevice.html#brls0yq
player = dsp.AudioPlayer('SampleRate',fs, 'DeviceDataType', dtyp); 
player.DeviceName = 'ASIO Fireface USB';

% Setting the output channels
channel_exc = 6;

% initialize parameter for excitation signal
rectime = 20; % recording time


amp1 =0.8;
% fss=2534;
% fss=7552;
% fss=7602;
% fss=2517;
fss=7422/3;   

fres = fss; %7422 2474
sample_time     = 1/fs;        
sample_frame    = fs*rectime/frameSize; % This is a sample per frame. 
% 

% Creation of Excitation Signal (SIN function)
sine1                   = dsp.SineWave;     %Create SinveWave Object
sine1.Frequency         = fres;
sine1.Amplitude         = amp1;  
sine1.PhaseOffset       = 0;
sine1.SampleRate        = fs;       
sine1.SamplesPerFrame = rectime/sample_time; % fs*rectime
%In the above line we dont devide the parameter to framezie so it means that 
%we generate one big sampe data on single frame (at DSP) and then we devide
%it to our desired framesize and send it to output channel.

amp2 =0.00001;
amp2 =amp1;

sine2                   = dsp.SineWave;     %Create SinveWave Object
sine2.Frequency         = fres;
sine2.Amplitude         = amp2;  
sine2.PhaseOffset       = 100;
sine2.SampleRate        = fs;       
sine2.SamplesPerFrame = rectime/sample_time;

% figure(1)
% plot(sine1());
% xlim([0 120])
% ylim([-1.2 1.2])
% title('Sin zero phase')
% figure(2)
% plot(sine2());
% xlim([0 120])
% ylim([-1.2 1.2])
% title('Sin with 179 degree phase shift')
% figure(3)
% plot(sine1()+sine2());
% xlim([0 120])
% ylim([-1.2 1.2])
% title('Sin sumation')



% Creation of Excitation Signal (Chirp function)
startfreq =7000;
endfreq = 8000;
endfreqedit= endfreq-startfreq;

chirpel                       = dsp.Chirp;    % Create Chirp Wave Object
chirpel.Type                  ='Linear'; % Specify the frequency sweep type as Swept cosine, Linear, Logarithmic, or Quadratic
chirpel.SweepDirection        ='Unidirectional'; %Specify the sweep direction as either Unidirectional or Bidirectional.
chirpel.TargetFrequency       =endfreq;
chirpel.InitialFrequency      =startfreq;
chirpel.OutputDataType        = 'double';
chirpel.TargetTime            =rectime;
chirpel.SweepTime             =rectime;
chirpel.SamplesPerFrame       =fs*rectime;
chirpel.SampleRate            =fs;

% Creation of Square Excitation Signal (pulse wave)

% fsp=96000;  % sampling rate
% endt=3; % generate a signal up to endt second
% t=0:1/(fsp):endt;   % Time span
% freq=2500; % frequency
% pulse1= 0.1*square(2*pi*freq*t);
% plot(t,x1);
% axis([0 5/freq -1.2 1.2]); % limiting the axis for better visualization


% Guasian pusle parameters:
fc = 7422; % Actual Frequency 
pd=0.3*rectime;% Control the Pulse duration
%-----------------------------------------
ts = (0:sample_time:rectime);  %time base
t0 = max(ts)/2; % Used to centering the pulse
gauss = (exp(-2*log(2)*(ts-t0).^2/(pd)^2)).*cos(-2*pi*fc*(ts-t0));
gauss=gauss';


% Select excitation type
% signal_exc = sine1();
% signal_exc = sine1()+sine2();
signal_exc = 1.05*chirpel();
%signal_exc = pulse1();
% signal_exc = 1.1*gauss();

% Allocating space for Recording Input/Output signals
sig = zeros(length(signal_exc),length(channel_exc));
exc = zeros(frameSize, 1);   %array, 'audioOut' to be played by the player is predefined

for start = 1:frameSize:length(signal_exc)-frameSize
    range = start:start+frameSize-1;
    
     % Play or send excitation signal
     step(player, exc);     
     exc(:,channel_exc)= signal_exc(range,1);
     %exc(:,channel_exc)= signal_exc(range);

     % Record   ... edit the following for recording the output in input 
     audioIn = step(recorder);                  % in 'audioIn' the current frame is saved for each pass
     sig(range,1)=audioIn(:,input_channel);    % 'in' stores array of measured data

% Y = fft(sig);
% L =fs*rectime;
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% plot(f,P1) 
% xlim([1000 8000])
% drawnow limitrate nocallbacks
end

%%
figure
plot (sig)
% ylim([-0.003 0.003])
% xlim([3.382*100000 3.392*100000])
title('Recorded Signal - fs = 96000, Rectime = 20, freq=7k-8k, Exc. 442mA')
% figure
% plot (signal_exc)
% ylim([-0.5 0.5])
% title('Excitation Signal')

%%
figure
%hold on
Y = fft(sig);
L =fs*rectime;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f,P1) 
xlim([1000 23000])
title('FFT responce of ME sensor')
xlabel('frequency (Hz)')
ylabel('Amplitude')
ydb = mag2db(P1);
figure
plot(f,ydb) 
xlim([1000 23000])
title('FFT responce of ME sensor')
xlabel('frequency (Hz)')
ylabel('Amplitude (db)')
% set(gca, 'YScale', 'log')
%%
%pause(har.QueueDuration);  % Wait until audio records to the end
release(player);
release(recorder);

