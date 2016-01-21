% official code to cut large 15mins hdf5 to 1 min hdf5 for every 7 min
% video
%
%
%
%
%

%% nas207-1
% from_pc207-7\ copied_from_pc207-13        copied_from_pc207-8
% from pc220-7
% from pc220-6
% from pc207-18

%% nas207-3
% pc207-8


%%

folder = 'pc207-14-Laura\';
path = ['N:\Kezhi\DataSet\AllFiles\MaskedVideos\nas207-3\',folder];

root_folder = genpath([path,'.']);

file=dir([path,'*.hdf5']);
num_file = size(file,1);

for nf = 1:num_file;
    
%hdf5_file = '247 JU438 on food L_2011_03_03__11_18___3___1';
%hdf5_file = '431 JU298 on food L_2011_03_03__12_29___3___5';
hdf5_file = file(nf).name(1:end-5);

hdf5_path = [path,hdf5_file,'.hdf5'];



% read key data
mask_info = h5info(hdf5_path, '/mask');
frame_size = mask_info.Dataspace.Size(1:2);
frame_total = mask_info.Dataspace.Size(3);
frame_pos = h5read(hdf5_path,'/vid_frame_pos'); % start from 0
time_pos = h5read(hdf5_path,'/vid_time_pos');

normalize_val = 1000;
% normalize_val = 1;
% while(1)
%    if time_pos(2)> 100
%        time_pos = time_pos/10;
%        normalize_val = normalize_val*10;
%    else
%        break
%    end
% end

video_length = time_pos(end)/normalize_val;
video_min = round(video_length/60);

if video_min > 7
    vv = round(video_min/15*2);
else
    vv = 1;
end
if vv> 10e5
    vv
    nf
    hdf5_file
    delete([path,hdf5_file,'.hdf5']);
    delete(['N:\Kezhi\DataSet\AllFiles\VideosExl\nas207-3\pc207-14-Laura\',hdf5_file,'.avi']);
    delete(['N:\Kezhi\DataSet\AllFiles\VideosExl\nas207-3\pc207-14-Laura\',hdf5_file,'.log.csv']);
    continue;
    %delete(['N:\Kezhi\DataSet\AllFiles\VideosExl\nas207-3\pc207-14-Laura',hdf5_file,'.hdf5']);
    %delete(['N:\Kezhi\DataSet\AllFiles\VideosExl\nas207-3\pc207-14-Laura',hdf5_file,'.log.csv']);    
   % error('vv is too large');
end
if video_min > 0
    randn_num = randperm(video_min);
else 
    randn_num = 1;
end
randn_start = randn_num(1:vv)-1;
start_time = randn_start*60;
recording_time = 60;   % 60

curr_root = [path,hdf5_file];
mkdir([curr_root,'_tif']);

for nn = 1:vv;
    
    sub_ind = randn_start(nn);
    % find logical indexes that we have interests, be careful '>=' and '<'
    logical_ind = (time_pos>= start_time(nn)*normalize_val & time_pos< (start_time(nn)+recording_time)*normalize_val);
    % get correspoinding time stamp
    time_vec = time_pos(logical_ind)/normalize_val;
    % find the corresponding frame index, +1 is because it starts from 1
    frame_ind = frame_pos(logical_ind) + 1;
    % when timestamp = 0, there is no corresponding frame
    if size(frame_ind)<1
        error('time stamp error');
    end
    if frame_ind(end) > frame_total
        frame_ind = frame_ind(1): frame_total;
        time_vec = time_vec(1:size(frame_ind,2));
    end
    
    %% create matrices that will be saved
    % build corresponding sub mask
    sub_mask = h5read(hdf5_path,'/mask',[1,1,double(frame_ind(1))],[frame_size(1),frame_size(2),double(frame_ind(end)-frame_ind(1)+1)]);
    substitute_intensity = round(median(median(sub_mask(sub_mask>0)))*1.1);
    sub_mask(sub_mask<1)=substitute_intensity;
    % save start time stamp
    time_start = min(time_vec);
    % save relative time stamp for frames
    sub_time = time_vec - time_start;
    
    
    %% save to hdf5
    output_file = [curr_root, '_tif','\',hdf5_file, '(',num2str(sub_ind),')','.hdf5'];
    
    h5create(output_file,'/mask',size(sub_mask),'Datatype','uint8','ChunkSize',[640 480 1],'Deflate', 4, 'Fletcher32',1,'Shuffle', 1);
    h5write(output_file, '/mask', sub_mask);
    h5create(output_file,'/vid_time_pos',size(sub_time),'FillValue',0);
    h5write(output_file, '/vid_time_pos', sub_time);
    h5create(output_file,'/time_start',1,'FillValue',0);
    h5write(output_file, '/time_start', time_start);
    
    % debug use
    if size(sub_mask,3)~=size(sub_time,1)
        sub_mask;
        sub_time
    end
    
    
    %% create tif files
    % write the first frame tif
    ii = 1;
    curr_img_name = [hdf5_file,'(',num2str(randn_start(nn)),').tif'];
    imwrite((sub_mask(:,:,ii))',[curr_root,'_tif','\',curr_img_name]);
    
    % write all other frames
    pre_timeStamp = 0;
    total_frame = size(sub_time,1);
    while(ii<total_frame)
        ii = ii + 1;
        cur_timeStamp = sub_time(ii);
        % when current time stamp increase to anther second, save the frame
        if floor(cur_timeStamp)>floor(pre_timeStamp);
            imwrite((sub_mask(:,:,ii))',[curr_root,'_tif','\',curr_img_name],'WriteMode','append');
        end
        pre_timeStamp = cur_timeStamp;
    end

end
end









