clear all
close all
%% params
rootdir = 'D:\Projects\SpikeSorting\data';
save_path='D:\Projects\SpikeSorting\figure\final';
min_frq=[1 ,4 ,30, 50];
max_frq=[4, 12, 50, 100];
x_tg=cell(0,numel(min_frq));
x_wt=cell(0,numel(min_frq));
y_tg=cell(0,numel(min_frq));
y_wt=cell(0,numel(min_frq));
tg_names={};
wt_names={};
%% select wt rats and channels
% removeing 3 top channels
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
% removeing 3 top channels
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
%% welech for wt
c=1; %counter
for num_wt=1:length(wt_files)
    files = dir(fullfile(rootdir, ['\**\final\*' wt_files(num_wt).day '*' wt_files(num_wt).name '*.mat']));
    for n=1:length(files)
        load([files(n).folder '\' files(n).name]);
        data = data';
        for channel=wt_files(num_wt).channels
            % saveing names for ploting
            wt_names{c}=replace(replace([files(n).name(1:end-4) ' ch ' num2str(channel)],'_',' '),'','');
            % compute pwlech
            [x_psd,y_psd]=single_channel_spectrum(data(channel,:));

            % extracting normalizition of 0.5-1 Hz
            y_n = y_psd(find(x_psd == 0.5) : find(x_psd == 1),1);
            y_psd=y_psd/mean(y_n);%normalizing


            for n_frq=1:numel(min_frq)
                x = x_psd(find(x_psd == min_frq(n_frq)) : find(x_psd == max_frq(n_frq)) ,1);
                y = y_psd(find(x_psd == min_frq(n_frq)) : find(x_psd == max_frq(n_frq)) ,1);
                x_wt(c,n_frq) = {x};
                y_wt(c,n_frq) = {y};
            end
            c=c+1;
        end
    end
end
%% welech for tg
c=1; %counter
for num_tg=1:length(tg_files)
    files = dir(fullfile(rootdir, ['\**\final\*' tg_files(num_tg).day '*' tg_files(num_tg).name '*.mat']));
    for n=1:length(files)
        load([files(n).folder '\' files(n).name]);
        data = data';
        for channel=tg_files(num_tg).channels
            % saveing names for ploting
            tg_names{c}=replace(replace([files(n).name(1:end-4) ' ch ' num2str(channel)],'_',' '),'','');
            % compute pwlech
            [x_psd,y_psd]=single_channel_spectrum(data(channel,:));

            % extracting normalizition of 0.5-1 Hz
            y_n = y_psd(find(x_psd == 0.5) : find(x_psd == 1),1);
            y_psd=y_psd/mean(y_n);%normalizing

            for n_frq=1:numel(min_frq)
                x = x_psd(find(x_psd == min_frq(n_frq)) : find(x_psd == max_frq(n_frq)) ,1);
                y = y_psd(find(x_psd == min_frq(n_frq)) : find(x_psd == max_frq(n_frq)) ,1);
                x_tg(c,n_frq) = {x};
                y_tg(c,n_frq) = {y};
            end
            c=c+1;
        end
    end
end
%% Ploting
pvalues=[];
for n_frq=1:numel(min_frq)
    close all
    y_tg_frq=cell2mat(y_tg(:,n_frq)');
    y_wt_frq=cell2mat(y_wt(:,n_frq)');

    X = categorical({'TG','WT'});
    X = reordercats(X,{'WT','TG'}); 
    b=bar(X(1),mean(mean(y_tg_frq)),'r');
    hold on

    bar(X(2),mean(mean(y_wt_frq)),'b');
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'fontweight','bold','fontsize',20);

    str = '#D9DDD6'; % colorrr
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    set(gcf,'Color',color); set(gca,'Color',color); set(gcf,'InvertHardCopy','off');
    
    %% saving
    if not(isfolder([save_path]))
        mkdir([save_path])
    end
    saveas(gcf,[save_path '\' num2str(min_frq(n_frq)) '_to_' num2str(max_frq(n_frq)) 'hz' '.svg']);
end
