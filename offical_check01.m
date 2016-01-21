% test if head/tail points are correct
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

root = 'C:\Kezhi\Videos\ReadyForLable\Andre\03-03-11\';
folder_name = '247 JU438 on food L_2011_03_03__11_18___3___1';

file=dir([root,folder_name,'\','*).tif']);
file_xls = dir([root,folder_name,'\','*).xlsx']);

num_file = size(file,1);
for nf = 1:2;
    %nf = 1;

    file_root = [root,folder_name,'\',file(nf).name];
    info = imfinfo(file_root);
    num_images = numel(info);

    xls_root = [root,folder_name,'\',file_xls(nf).name];
    current_xls = xlsread(xls_root);
    if max(current_xls(:,6))>640 || max(current_xls(:,7))>480
        current_xls(:,6) = current_xls(:,6)/225777.78*640;
        current_xls(:,7) = current_xls(:,7)/169333.33*480;
    end

     for k = 1:num_images;
        A = imread(file_root, k);
        mm = k *2 -1;
       figure(10), imshow(A);
       hold on
       plot(current_xls(mm,6),current_xls(mm,7),'r*');
       plot(current_xls(mm+1,6),current_xls(mm+1,7),'r*');
       hold off
       pause(0.5);
     end
end
 