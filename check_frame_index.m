clear
clc

path = 'C:\Kezhi\MyCode!!!\ManualVideos\';

% please add the folder name here
addpath(genpath([path,'.']));

ffmpeg = 'C:\FFMPEG\bin\ffmpeg';
ffprobe = 'C:\FFMPEG\bin\ffprobe';

name = '135 CB4852 on food L_2011_03_09__15_51_36___1___8';
folder = '09-03-11\';
root = ['N:\Kezhi\DataSet\ReadyForLable\nas207-1\from_pc207-7\copied_from_pc207-13\Andre\',folder];
subfolder = [root, name, '\'];

video_full = videoread([root,'135 CB4852 on food L_2011_03_09__15_51_36___1___8.avi']);
video_mm = mmread([root,'135 CB4852 on food L_2011_03_09__15_51_36___1___8.avi']);
video_rV = videoReader([root,'135 CB4852 on food L_2011_03_09__15_51_36___1___8.avi']);
video_RV = VideoReader([root,'135 CB4852 on food L_2011_03_09__15_51_36___1___8.avi']);


vv = 15;
curr_ind = 1;
for nn = 1:vv;
    curr_file = [root,name,'_avi', '\',name, '(',num2str(nn-1),')','.avi'];
    curr_video = mmread(curr_file);   %curr_video = videoread(curr_file);
    
    num_ind = size(curr_video.frames,2);
    
    frame_diff = sum(sum(sum( abs(curr_video.frames(1).cdata - video_full.frames(curr_ind).cdata))))
    curr_ind
    curr_ind = num_ind + curr_ind;
    curr_ind-1
end