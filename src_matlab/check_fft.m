clear;
close all;
format long;
clc;
Fs = 10; % sampling frequency 1 kHz
t =  [0,10,20,30,40,50,60,70,80,90]; % time scale
x = [10,120,130,120,120,100,123,456,78,89]; % time series
x=x-mean(x);
plot(t,x), axis('tight'), grid('on'), title('Time series'), figure
nfft = 512; % next larger power of 2
y = fft(x,nfft); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
%y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
%f_scale = (0:nfft/2)* Fs/nfft; % frequency scale
f_scale = (0:nfft-1)*Fs/nfft;
plot(f_scale, y),axis('tight'),grid('on'),title('Dominant Frequency')
fest = f_scale(k); % dominant frequency estimate
fprintf('Dominant freq.: true %f Hz, estimated %f Hznn', f, fest)
fprintf('Frequency step (resolution) = %f Hznn', f_scale(2))