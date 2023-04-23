rootdir = 'D:\Projects\SpikeSorting\data';
save_path='D:\Projects\SpikeSorting\report\bz_coupling\';
filelist = dir(fullfile(rootdir, '\**\final\**\*.mat'));
for n=1:length(filelist)
    clearvars -except rootdir save_path filelist n
    close all
    load([filelist(n).folder '\' filelist(n).name])
    samplingRate = 1000;            % Sampling frequency                    
    T = 1/samplingRate;             % Sampling period       
    N = height(data);             % Length of signal
    timestamps = (0:N-1)*T;        % Time vector
    timestamps=timestamps';
    lfp.data=data;
    lfp.timestamps=timestamps;
    lfp.samplingRate=samplingRate;
    lfp.channels=[1:14]';
    [comod] = bz_CrossFreqMod(lfp,[4:0.5:12],[1:1:140]);

        %% saving
    if not(isfolder([save_path]))
        mkdir([save_path])
    end
    saveas(gcf,[save_path filelist(n).name '.fig']);
end

