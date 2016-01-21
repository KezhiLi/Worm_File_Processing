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

% C:\Kezhi\Videos\ReadyForLable\Andre\
% 03-03-11
% 04-03-11
% 07-03-11
% 08-03-11
% 09-03-11
% 10-03-11
% 12-04-11
% 13-04-11
% 14-04-11
% 15-03-11
% 17-02-11
% 17-03-11

% C:\Kezhi\Videos\ReadyForLable\pc207-8-Laura

root = 'C:\Kezhi\Videos\ReadyForLable\pc207-8-Laura\03-08-10\';

file=dir([root,'*.avi']);
% for n=1:length(file)
%     temp=dlmread(['E:\new\',file(n).name],' ',0,1);
%     eval([file(n).name(1:end-4),'=temp;'])
% end
num_file = size(file,1);

%for nf = 2:num_file;
nf = 9;

   %  video = mmread([root,file(nf).name],[],[0 900]); %get only the first 10 frames
    video = mmread([root,file(nf).name]); %get only the first 10 frames
    
    mkdir(root,file(nf).name(1:end-4));
    curr_root = [root,file(nf).name(1:end-4)];
    for ii = 1:size(video.frames,2);
        jj = floor((ii-1)/1800);
        kk = mod(ii, 1800);
        if kk == 1;
            curr_img_name = [file(nf).name(1:end-4),'(',num2str(jj+1),').tif'];
            img = rgb2gray(video.frames(ii).cdata);
            imwrite(img,[curr_root,'\',curr_img_name]);
        else
            subsamp_ind = mod(kk,30);
            if  subsamp_ind == 1
                img = rgb2gray(video.frames(ii).cdata);
                imwrite(img,[curr_root,'\',curr_img_name],'WriteMode','append');
            end
        end
    end
    
    clearvars video
%end



%movie(video.frames)
%
% nFrames = 300;
% mov(1:nFrames) = struct('cdata', [],'colormap', []);
% for k = 1:nFrames
%    imshow(video.frames(k).cdata);
%    mov(k) = getframe(gcf);
% end
%
% % Create AVI file.
% movie2avi(mov, 'test1.avi', 'compression', 'None');