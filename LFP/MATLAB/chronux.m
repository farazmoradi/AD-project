clear all
%% params
rootdir = 'D:\Projects\SpikeSorting\data';
srate=1000;
min_frq=[1,4,30,50]; 
max_frq=[4,12,50,100];
clim_max=[3,8.5,0.4,0.25,8.5];
clim_min=[0,0,0,0,0];
numfrex=40;
pnts=900;
tf_avg = zeros(numfrex,pnts*srate);
save_path ='D:\Projects\SpikeSorting\figure\final\chronux\';
%% select wt rats and channels
wt_files(1).name='WT4';
wt_files(1).channels=[1,2,3,4,5,6,8,10,11,12,13];
wt_files(1).day='d1';
wt_files(2).name='WT4';
wt_files(2).channels=[1,2,3,4,5,6,8,10,11,12,13];
wt_files(2).day='d2';
wt_files(3).name='WT4';
wt_files(3).channels=[1,2,3,4,5,6,8,10,11,12,13];
wt_files(3).day='d3';
wt_files(4).name='WT5';
wt_files(4).channels=[1,2,3,4,5,6,8,10,11,12,13];
wt_files(4).day='d1';
wt_files(5).name='WT5';
wt_files(5).channels=[1,2,3,4,5,6,8,10,11,12,13];
wt_files(5).day='d2';
wt_files(6).name='WT6';
wt_files(6).channels=[1,2,3,4,5,6,8,10,11,12,13];
wt_files(6).day='d1';
wt_files(7).name='WT6';
wt_files(7).channels=[1,2,3,4,5,6,8,10,11,12,13];
wt_files(7).day='d2';
%% selecting tg rats and channels
tg_files(1).name='TG1';
tg_files(1).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(1).day='d1';
tg_files(2).name='TG2';
tg_files(2).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(2).day='d1';
tg_files(3).name='TG2';
tg_files(3).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(3).day='d2';
tg_files(4).name='TG2';
tg_files(4).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(4).day='d3';
tg_files(5).name='TG3';
tg_files(5).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(5).day='d1';
tg_files(6).name='TG3';
tg_files(6).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(6).day='d2';
tg_files(7).name='TG3';
tg_files(7).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(7).day='d3';
tg_files(8).name='TG4';
tg_files(8).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(8).day='d1';
tg_files(9).name='TG4';
tg_files(9).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(9).day='d2';
tg_files(10).name='TG5';
tg_files(10).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(10).day='d1';
tg_files(11).name='TG5';
tg_files(11).channels=[1,2,3,4,5,6,8,10,11,12,13];
tg_files(11).day='d2';
%% run wavelet
counter_avg=0;
params.tapers=[5 9];
params.Fs = 1e3;
params.trialave =0;
movingwin =[15 15];   %  in the form [window winstep] 
                       %  i.e length of moving window and step size
if not(isfolder([save_path]))
    mkdir([save_path])
end
%% time-frq wt
for c=1:numel(min_frq)
    close all
    S1 =[];
    for num_wt=1:length(wt_files)
        sub_files = dir(fullfile(rootdir, ['\**\final\*' wt_files(num_wt).day '*' wt_files(num_wt).name '*.mat']));
        for n=1:length(sub_files)
            load([sub_files(n).folder, '\' ,sub_files(n).name])
            for channel=wt_files(num_wt).channels
                params.fpass= [min_frq(c) max_frq(c)];
                signal = data(1:900*srate,channel)';
                [S,t,f] = mtspecgramc( signal, movingwin, params );
                params.fpass= [0.5 1];
                [S_mean,temp1,temp2] = mtspecgramc( signal, movingwin, params );
                baseline=mean(S_mean,2);
                S=S./baseline;
                S1 = cat(3,S1,S); 
            end
            clear data
        end
    end
    close all
    S_final= mean(S1,3);
    %% ploting
    contourf(t,f,S_final',40,'linecolor','none')
    set(gca,'fontweight','bold','fontsize',40)
    set(gcf, 'Position', get(0, 'Screensize'));
    colormap(jet)
    colorbar
    set(gca,'clim',[clim_min(c) clim_max(c)]);

    str = '#D9DDD6'; % color
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    set(gcf,'Color',color); set(gca,'Color',color); set(gcf,'InvertHardCopy','off');
    saveas(gcf,[save_path num2str(min_frq(c)) 'Hz to ' num2str(max_frq(c)) 'Hz_wt.svg']);
end
%% time-frq tg
for c=1:numel(min_frq)
    close all
    S1 =[];
    for num_tg=1:length(tg_files)
        sub_files = dir(fullfile(rootdir, ['\**\final\*' tg_files(num_tg).day '*' tg_files(num_tg).name '*.mat']));
        for n=1:length(sub_files)
            load([sub_files(n).folder, '\' ,sub_files(n).name])
            for channel=tg_files(num_tg).channels
                params.fpass= [min_frq(c) max_frq(c)];
                signal = data(1:900*srate,channel)';
                [S,t,f] = mtspecgramc( signal, movingwin, params );
                params.fpass= [0.5 1];
                [S_mean,temp1,temp2] = mtspecgramc( signal, movingwin, params );
                baseline=mean(S_mean,2);
                S=S./baseline;
                S1 = cat(3,S1,S); 
            end
            clear data
        end
    end
    close all
    S_final= mean(S1,3);
    %% ploting
    contourf(t,f,S_final',40,'linecolor','none')
    set(gca,'fontweight','bold','fontsize',40)
    set(gcf, 'Position', get(0, 'Screensize'));
    colormap(jet)
    colorbar
    set(gca,'clim',[clim_min(c) clim_max(c)]);
    str = '#D9DDD6'; % color
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    set(gcf,'Color',color); set(gca,'Color',color); set(gcf,'InvertHardCopy','off');
    saveas(gcf,[save_path num2str(min_frq(c)) 'Hz to ' num2str(max_frq(c)) 'Hz_tg.svg']);
end
