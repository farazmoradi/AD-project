clear all
rootdir = 'D:\Projects\SpikeSorting\data';
save_path='D:\Projects\SpikeSorting\data\mi_index\1_17_2023\';
filelist = dir(fullfile(rootdir, '\**\final\**\*.mat'));
mi_table=table;
for n=1:length(filelist)
    clearvars -except rootdir save_path filelist n mi_table
    close all
    load([filelist(n).folder '\' filelist(n).name])
    srate = 1000;            % Sampling frequency
    mi_mat=[];
    p_mat=[];
    data(:,[9,7,14])=[];
    data=data';
    subject_mi=[];
    subject_amp=[];
    subject_mraw=[];
    for channel=1:11
        close all
        clearvars -except subject_mraw subject_mi subject_amp srate data channel n filelist save_path rootdir mi_mat p_mat mi_table
        
        dt = 1/srate;
        t = dt*(1:length(data(channel,:)));
        
        x = data(channel,:);
        
        %%
        numpoints=length(x); %% number of sample points in raw signal
        numsurrogate=200; %% number of surrogate values to compare to actual value
        minskip=srate; %% time lag must be at least this big
        maxskip=numpoints-srate; %% time lag must be smaller than this
        skip=ceil(numpoints.*rand(numsurrogate*2,1));
        skip(find(skip>maxskip))=[];
        skip(find(skip<minskip))=[];
        skip=skip(1:numsurrogate,1);
        surrogate_m=zeros(numsurrogate,1);
        %% HG analytic amplitude time series, uses eegfilt.m from EEGLAB toolbox
        %% (http://www.sccn.ucsd.edu/eeglab/)
        amplitude=abs(hilbert(eegfilt(x,srate,50,100)));
        %% theta analytic phase time series, uses EEGLAB toolbox
        phase=angle(hilbert(eegfilt(x,srate,4,12)));
        %% complex-valued composite signal
        z=amplitude.*exp(i*phase);
        %% mean of z over time, prenormalized value
        m_raw=mean(z);
        %% compute surrogate values
         for s=1:numsurrogate
         surrogate_amplitude=[amplitude(skip(s):end) amplitude(1:skip(s)-1)];
         surrogate_m(s)=abs(mean(surrogate_amplitude.*exp(i*phase)));
         disp(surrogate_m(s))
         disp(numsurrogate-s)
         end
        %% fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
        [surrogate_mean,surrogate_std]=normfit(surrogate_m);
        %% normalize length using surrogate data (z-score)
        m_norm_length=(abs(m_raw)-surrogate_mean)/surrogate_std;
        m_norm_phase=angle(m_raw);
        m_norm=m_norm_length*exp(i*m_norm_phase);
%% saving
        t_subject=cell2table({[filelist(n).name 'channel' num2str(channel)],abs(m_norm),abs(m_raw),mean(amplitude)},'VariableNames',{'name','mi_index_norm','m_raw','gamma'})
        mi_table=[mi_table;t_subject]
    end

%     saveas(gcf,[save_path filelist(n).name '.fig']);
end
    %% saving
if not(isfolder([save_path]))
    mkdir([save_path])
end
writetable(mi_table,[save_path 'mi_table_50_100_revised.csv'],'Delimiter',',','QuoteStrings',true)

% figure(3)
% bar(10:20:710,[table2array(mi_table(17,3)) table2array(mi_table(17,3))])
% xlabel('\Theta Phase (Deg)')
% ylabel('Normalized mean \gamma Amplitude')
% ylim([0 0.1])
% set(gca,'xtick',0:90:720)
