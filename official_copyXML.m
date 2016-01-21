% script to copy all .info.xml from Z:\ to N:\, in corresponding folder.
% Change path and xml_path to the folder that is going to be copied
%
% 
% 
% the task has completed on 02/12/2015
% Kezhi Li, Imperial College London

folder = 'copied_from_pc207-13\';
path = ['N:\Kezhi\DataSet\AllFiles\VideosExl\nas207-1\from pc207-18\'];
xml_path = ['Z:\thecus\nas207-1\experimentBackup\from pc207-18\!worm_videos\'];

root_folder = genpath([path,'.']);

% find all .avi files
file=dir([path,'*.avi']);
% number of avi files in path folder
num_file = size(file,1);

for nf = 1:num_file;
    % show process
    nf
    % file name of avi
    avi_name = file(nf).name(1:end-4);

    % xml file name
    xml_name = [avi_name,'.info.xml'];
    %[status,list]=system(['dir ',xml_path,'*',xml_name]);
    % find the xml file in all folder and subfolder
    xml_file = subdir([xml_path, '*',xml_name]);
    
    % see if corresponding file exists
    if size(xml_file,1) ==1 
        copyfile(xml_file.name,path);
    elseif size(xml_file,1) >1 
        copyfile(xml_file(1).name,path);
    else
        msg = ['can not find file',xml_name]
    end
    
end