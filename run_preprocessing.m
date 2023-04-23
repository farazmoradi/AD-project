clear all
clc
rootdir = 'D:\Projects\SpikeSorting\data';
filelist = dir(fullfile(rootdir, '\**\*.edf')); 
for n=1:length(filelist)
   preprocessing(filelist(n).folder,filelist(n).name)
end
