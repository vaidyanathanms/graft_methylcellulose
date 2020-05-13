function [outfft,freq_out] = compute_fft(indata, yindim, tsample, nsamples)
% Store N-equispaced samples from indata in the interval tminval and tmaxval
% and compute fft and return
% ref page: https://www.mathworks.com/matlabcentral/answers/155036-how-to-plot-fft-of-time-domain-data
% ref page: https://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/
%% Load data
% store sampled values in a separate file

invals = indata(:,yindim);  % Load kappasq data
Ts = tsample;               % Sampling Interval (tau)
Fs = 1/Ts;                  % Sampling Frequency (tau^-1)
invals(isnan(invals))=[];   % Eliminate ‘NaN’ Values First
len_inval = size(invals,1); % Length of Data Vector

%% Fourier Transform

invals   = invals - mean(invals);                % subtract DC component
outfft   = 1/len_inval*fft(invals,nsamples);     % Fourier Series of input data: nsamples is the # of sampling points for FFT
freq_out = linspace(0,1,fix(len_inval/2)+1)'*Fs; % output one sided frequency vector
outfft(2:end) = 2*outfft(2:end);                 % except first value (DC) all other values are multiplied by 2;

end

