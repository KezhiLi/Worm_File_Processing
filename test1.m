% 
% 
% 
% 
% 
% 
clear
clc

path = 'C:\Kezhi\MyCode!!!\ManualVideos\';

% please add the folder name here
addpath(genpath([path,'.']));


root = 'C:\Kezhi\MyCode!!!\Sample\video\Real_data\experimentBackupfrom pc207-7worm_videoscopied from pc207-8Andre03-03-11\';

file=dir('E:\new\*.txt');
for n=1:length(file)
    temp=dlmread(['E:\new\',file(n).name],' ',0,1);
    eval([file(n).name(1:end-4),'=temp;'])
end


video = mmread([root,'247 JU438 on food L_2011_03_03__11_18___3___1.avi'],18901:19200); %get only the first 10 frames

%movie(video.frames)

nFrames = 300;
mov(1:nFrames) = struct('cdata', [],'colormap', []);
for k = 1:nFrames 
   imshow(video.frames(k).cdata);
   mov(k) = getframe(gcf);
end

% Create AVI file.
movie2avi(mov, 'test1.avi', 'compression', 'None');