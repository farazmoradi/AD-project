rootdir = 'D:\Projects\SpikeSorting\data';
save_path='D:\Projects\SpikeSorting\report\hvs\';
filelist = dir(fullfile(rootdir, '\**\final\**\*.mat'));
hvs_table=table;
for n=1:length(filelist)
    clearvars -except rootdir save_path filelist n hvs_table
    load([filelist(n).folder '\' filelist(n).name])
    data(:,[9,7,14])=[];
    srate=1000;
    % bandpass 6 12hz
    order=3;
    fl=6;
    fh=12;
    wn= [fl fh] / (srate/2);
    type= 'bandpass';
    [b,a]= butter(order,wn,type);
    data_6_12= filtfilt(b,a,data);
    
    % bandpass 1 4hz
    order=3;
    fl=1;
    fh=4;
    wn= [fl fh] / (srate/2);
    type= 'bandpass';
    [b,a]= butter(order,wn,type);
    data_1_4= filtfilt(b,a,data);
    
    h_data_6_12=mean(abs(hilbert(data_6_12)).^2,2);
    h_data_1_4=mean(abs(data_1_4).^2,2);
    
    ratio=h_data_6_12./h_data_1_4;
    greater_than_2=ratio>2;
    theta_epoch=[];
    counter=1;
    i=1;
    while i<=length(greater_than_2)-2000
        s=0;
        while greater_than_2(i)==1 && i<=length(greater_than_2)-2000
            s=s+1;
            i=i+1;
        end
        i=i+1;
        if s>=2000
            theta_epoch=[theta_epoch,{i-s-1;s-1}];
            counter=counter+1;
        end
    end
    theta_epoch=cell2mat(theta_epoch);
    
    avg_data = mean(data,2);
    idx=[]
    for i=1:length(theta_epoch)
        theta=avg_data(theta_epoch(1,i):theta_epoch(1,i)+theta_epoch(2,i));
        theta=theta>0.3;
        sum(theta)
        if sum(theta)>=1
            idx=[idx;i];
        end
    end
    hvs=length(idx)/(length(avg_data)/(srate*60));
    t_subject=cell2table({filelist(n).name,hvs},'VariableNames',{'name','hvs'})
    hvs_table=[hvs_table;t_subject]
end
    %% saving
if not(isfolder([save_path]))
    mkdir([save_path])
end
writetable(hvs_table,[save_path 'hvs.csv'],'Delimiter',',','QuoteStrings',true)
