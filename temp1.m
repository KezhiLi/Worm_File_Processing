path = 'C:\Kezhi\MyCode!!!\ManualVideos\';

% please add the folder name here
addpath(genpath([path,'.']));



ffprobe = 'C:\FFMPEG\bin\ffprobe';

folder = 'N:\Kezhi\DataSet\ReadyForLable\nas207-1\from_pc207-7\copied_from_pc207-8\Andre\03-03-11\431 JU298 on food L_2011_03_03__12_29___3___5_avi\';
name = '431 JU298 on food L_2011_03_03__12_29___3___5(9)';
input_file = [folder,name, '.avi'];
input_file_com = ['"' input_file '"'];
cmd_info = sprintf('%s -i %s -show_entries format=duration -v quiet -of csv="p=0"', ffprobe, input_file_com);
[status,video_length] = system(cmd_info)



 %video = mmread(input_file);
 video_full = mmread('N:\Kezhi\DataSet\ReadyForLable\nas207-1\from_pc207-7\copied_from_pc207-8\Andre\03-03-11\431 JU298 on food L_2011_03_03__12_29___3___5.avi');
 
 frame_file = 1801.8;
 %sum(sum(rgb2gray(video.frames(1).cdata)- rgb2gray(video_full.frames(frame_file*2+1).cdata)))
 
curr_frame = 0;
prev_num_frame = 0; 
for ii = 0:14;
name = ['431 JU298 on food L_2011_03_03__12_29___3___5(',num2str(ii),')'];
input_file = [folder,name, '.avi'];
input_file_com = ['"' input_file '"'];
video = mmread(input_file);

curr_frame = curr_frame + prev_num_frame;

sum(sum(rgb2gray(video.frames(1).cdata)- rgb2gray(video_full.frames(round(curr_frame+1)).cdata)))
result =[curr_frame+1,curr_frame+size(video.frames,2)]
prev_num_frame = size(video.frames,2);
end
