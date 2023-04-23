clear all
close all
%% params
rootdir = 'D:\Projects\SpikeSorting\data';
save_path=['D:\Projects\SpikeSorting\report\'];
min_frq=[1,1 ,4 ,30, 50];
max_frq=[100,4, 12, 50, 100];
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
%% plotting
for n_frq=1:numel(min_frq)
    close all
    x_tg_frq=cell2mat(x_tg(:,n_frq)');
    y_tg_frq=cell2mat(y_tg(:,n_frq)');

    x_wt_frq=cell2mat(x_wt(:,n_frq)');
    y_wt_frq=cell2mat(y_wt(:,n_frq)');
    % TG figure
    figure_tg=figure(1);
    plot(x_tg_frq,y_tg_frq)
    xlabel("Frequncy (Hz)")
    ylabel("Power")
    set(gca,'fontweight','bold','fontsize',16);
    legend(tg_names,'Location', 'Northeastoutside','FontSize',7)
    set(gcf, 'Position', get(0, 'Screensize'));
    xlim([min_frq(n_frq),max_frq(n_frq)])
    
    % WT figure
    figure_wt=figure(2);
    plot(x_wt_frq,y_wt_frq)
    xlabel("Frequncy (Hz)")
    ylabel("Power")
    set(gca,'fontweight','bold','fontsize',16);
    legend(wt_names,'Location', 'Northeastoutside','FontSize',8)
    set(gcf, 'Position', get(0, 'Screensize'));
    xlim([min_frq(n_frq),max_frq(n_frq)])

    
    %TG and WT figure
    figure_tg_wt=figure(3)
    
    avg_x_tg=mean(x_tg_frq');
    avg_y_tg=mean(y_tg_frq');
    
    avg_x_wt=mean(x_wt_frq');
    avg_y_wt=mean(y_wt_frq');
    
    std_tg=std(y_tg_frq');
    std_wt=std(y_wt_frq');
    
    %plotting tg
    
    curve1 = avg_y_tg + std_tg;
    curve2 = avg_y_tg - std_tg;
    x2 = [avg_x_tg, fliplr(avg_x_tg)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2', inBetween', 'r','FaceAlpha',0.1,'LineStyle','none');
    hold on;
    plot(avg_x_tg, avg_y_tg, 'r', 'LineWidth', 8);
    hold on;
    
    %plotting wt
    curve1 = avg_y_wt + std_wt;
    curve2 = avg_y_wt - std_tg;
    x2 = [avg_x_wt, fliplr(avg_x_wt)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'b','FaceAlpha',0.1,'LineStyle','none');
    hold on;
    plot(x_wt_frq, avg_y_wt, 'b', 'LineWidth', 8);
    legend("","TG","","WT",'color','none')
    xlabel("Frequncy (Hz)")
    ylabel("Power")
    set(gca,'fontweight','bold','fontsize',16);
    set(gcf, 'Position', get(0, 'Screensize'));
    xlim([min_frq(n_frq),max_frq(n_frq)])
    str = '#D9DDD6'; % color
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    set(gcf,'Color',color); set(gca,'Color',color); set(gcf,'InvertHardCopy','off');
    
    if not(isfolder([save_path]))
        mkdir([save_path])
    end
    
    saveas(figure_tg_wt,[save_path '\PSD_WT_TG_normed_0.5_1hz' num2str(min_frq(n_frq)) '_to_' num2str(max_frq(n_frq)) 'hz' '.svg']);
end