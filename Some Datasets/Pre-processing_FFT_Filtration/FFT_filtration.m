% Load the data_set for FFT 
load('10filtered_700.mat');
sig = v2' ;
%sig = datasel(:, 2);
%sig_new = reshape(sig(2:end), 1920000, 1);
% For Laser data

% Defining the parameters for FFT %
rect = 10;              %recorded time
Fs = 96000;           % Sampling frequency                    
T = 1 / Fs;             % Sampling period       
L = Fs * rect;          % Length of signal
%L = 1500;
%t = (0:L-1)* T;         % Time vector
%t = distancedataset(:, 1);
%t = (0:L)* T;

% For Gaussian pulse

% To improve the performance of fft, identify an input length
% that is the next power of 2 from the original signal length.
% Calling fft with this input length,
% pads the pulse X with trailing zeros to the specified transform length.
% n = 2^nextpow2(L);

% Signal filtration based on frequency (highpass,low pass, bandpass filters) %

high_filtered = highpass(sig(:,1),6000,Fs);
low_filtered  = lowpass(high_filtered, 38000, Fs);
% Y2 = bandpass(sig,[7422,14844],Fs);
% Y3 = bandpass(Y2,[14844,22266],Fs);
% Y4 = bandpass(sig,[22266,29688],Fs);
% Y5 = bandpass(sig,[29688,37110],Fs);


%Perform FFT %
Y1 = fft(sig); 
%Y1 = fft(low_filtered);
% Y1 = fft(displacement);
% Y1 = fft(simulated(:,1));

 % for gaussian pulse

%Y1 = fft(sig(:, 1), n);  
%Y1 = fft(low_filtered, n); 

% shift the zero frequency to the center %
 
%fshift = (-L/2:L/2)*(Fs/L); 
% fshift = (-L/2:L/2-1)*(Fs/L);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L. %
P2 = abs(Y1/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define Frequency domain. On average, longer signals produce better frequency approximations. %

f = Fs*(0:(L/2))/L;

% f = Fs*(0:(n/2))/n;  % for gaussian pulse
% P = abs(Y1/n).^2;

% plot(f,P(1:n/2+1)) 
% title("Gaussian Pulse in Frequency Domain")
% xlabel("f (Hz)")
% ylabel("|P(f)|^2")

figure;
plot(f,P1)
xlim([1000 40000])
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

hold on;      % hold the currenmt plot and dont overwrite

% The relationship between magnitude and decibels is ydb = 20 log10(y). %
figure;
ydb = mag2db(P1);
plot(f,ydb)
xlim([1000 40000])
title('FFT responce of velocity readout')
xlabel('frequency (Hz)')
ylabel('Amplitude (db)')

hold on;

% Finding peaks automatically. This will return the peak amplitudes in peaks and their corresponding frequencies in locs %
% 
% filtered_signal= abs(filtered_signal);
% [peaks, locs] = findpeaks(filtered_signal, fshift);
% 
% figure;
% % stem(locs, peaks, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r');
% stem(locs, peaks, 'filled');
% xlim([10 20000])
% title("Peaks and their locations  (frequencies)")
% xlabel('frequency (Hz)')
% ylabel('Amplitude')
% 
% hold off;

% Filtering based on the amplitudes of the frequencies %

% signal_fft_mag = abs(Y1);
% 
% threshold = power(4, -4);
% indices_to_filter = find(signal_fft_mag < threshold);
% Y1(indices_to_filter) = 0;
% 
% signal_filtered = ifft(Y1);
% 
% figure;
% 
% subplot(2,1,1);
% %sig = sig(10000:11000, 1);
% plot(t, sig);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Original Signal');
% 
% subplot(2,1,2);
% plot(t, signal_filtered);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Filtered Signal');
