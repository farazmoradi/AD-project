function preprocessing(folder,file_name)
%% Loading data
data = edfread([folder '\' file_name]);
data = cell2mat(data{:,:});
srate = 10000;            % Sampling frequency                    
T = 1/srate;             % Sampling period       
N = width(data);             % Length of signal
t = (0:N-1)*T;        % Time vector
hz = linspace(0,srate/2,floor(N/2)+1);
%% remove EMG and WIFI channels
data(:,17)=[];
data(:,11)=[];
data(:,5)=[];
%% remove fisrt and last 60s for possible artifacts
data=data(60*srate:end-60*srate,:);
%% bandpass filter
order=3;
fl=0.1;
fh=300;
wn= [fl fh] / (srate/2);
type= 'bandpass';
[b,a]= butter(order,wn,type);
data= filtfilt(b,a,data);
%% Downsampling
data = downsample(data,10);
%% Saving mat file
if not(isfolder([folder]))
    mkdir([folder])
end
save([folder file_name(1:end-4) '.mat'],'data')