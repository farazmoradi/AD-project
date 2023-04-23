function [x1,y1]=single_channel_spectrum(signal)
% signal: 1*time point
%% params
srate = 1000;
nfft = srate*100;
winsize = 2*srate; % 2-second window
hannw = .5 - cos(2*pi*linspace(0,1,winsize))./2;
%% pwelch
[y1,x1] = pwelch(signal,hannw,round(winsize/4),nfft,srate);
end
