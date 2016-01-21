clear
clc
% 
% nf = 7; % 1 
% 
% path2 = 'C:\Kezhi\MyCode!!!\ManualVideos\Check_Align_samples\';
% 
% % please add the folder name here
% addpath(genpath([path2,'.']));
% 
% ffmpeg = 'C:\FFMPEG\bin\ffmpeg';
% ffprobe = 'C:\FFMPEG\bin\ffprobe';
% 
% 
% root = 'C:\Kezhi\MyCode!!!\ManualVideos\Check_Align_samples\';
% folder = root;
% hdf_folder = folder;
% 
% % folder = 'MissingFrames_example\';
% % root = ['N:\Kezhi\DataSet\',folder];
% 
% root_folder = genpath([root,'.']);
% 
% file=dir([root,'*.hdf5']);
% num_file = size(file,1);
% 
% % for nf = 1:1;
% % choose index of file
% % nf = 3; % 1 
%     name  = file(nf).name(1:end-5);
%     excel_name = [name, '.log.csv'];
%     hdf5_name = [name, '.hdf5'];
% 
% % import csv data
% ss = importdata([folder,excel_name]);
% 
% % read data: real time, media time, and stage coordinates xy
% real_time = ss.textdata(2:end,1);
% media_time = ss.textdata(2:end,2);
% stage_xy = ss.data;
% 
% % check if number of rows in csv is equal to number of frames 
% data_rows = size(ss.data,1);
% if data_rows+1 ~= size(ss.textdata,1);
%     error('excel file problem: number of rows are not in uniform')
% end
% 
% % convert time from text to number 
% media_time_vec = [];
% 
% % convert time to seconds value
% for ii = 1:data_rows;
%     str1 = media_time(ii);
%     t1 = datevec(str1);
%     media_time_vec(ii) = t1(5)*60 + t1(6);
% end
% 
% hdf5_path = [hdf_folder,hdf5_name];
% 
% %% read key data
% % mask information
% mask_info = h5info(hdf5_path, '/mask');
% % mask matrix, 3 dimentions
% mask = h5read(hdf5_path, '/mask');
% % size of each frame
% frame_size = mask_info.Dataspace.Size(1:2);
% frame_total = mask_info.Dataspace.Size(3);
% frame_pos = h5read(hdf5_path,'/vid_frame_pos'); % start from 0
% time_pos = h5read(hdf5_path,'/vid_time_pos');
% normalize_val = 1000;
% if frame_total > size(time_pos)
%     error('number of frames is larger than number of time stamp');
% end
% 
% 
% %% calculate the central moving difference between each frame, and save them in 'diff_mask_central(:,1,2)'
% mask_bw = []; 
% CC_max = logical(zeros(frame_size));
% mask_central = zeros(frame_total,3);
% mask_previous = mask(:,:,1);
% 
% for ii = 1:frame_total;
%     ii
%     mask_backadjust = mask(:,:,ii);
%     substitute_intensity = round(median(median(mask_backadjust(mask_backadjust>0)))*1.1);
%     mask_backadjust(mask_backadjust<1)=substitute_intensity;
%     level = graythresh(mask_backadjust);
%     mask_bw = im2bw(mask_backadjust,level)<0.5;
%     CC = bwconncomp(mask_bw);
%     
%     numPixels = cellfun(@numel,CC.PixelIdxList);
%     [biggest,idx] = max(numPixels);
%     CC_max_this = CC_max;
%     CC_max_this(CC.PixelIdxList{1,idx}) = logical(1);
%     S = regionprops(CC_max_this,'Centroid');
%     mask_central(ii,1:2)=S.Centroid;
% 
%     if ii == frame_total
%         break;
%     else
%        pos_diff_pixel = double(mask(:,:,ii+1)) - double(mask(:,:,ii));
%         sum_postive = sum(sum(pos_diff_pixel(pos_diff_pixel>20)));
%         sum_negtive = sum(sum(pos_diff_pixel(pos_diff_pixel<-20)));
%        pixel_diff(ii,1:2) = [sum_postive, sum_negtive]; 
%     end
% end
% % normalize time
% mask_central(:,3) = time_pos(1:size(mask_central,1))/normalize_val;
% pixel_diff(:,3) = abs(pixel_diff(:,1).*pixel_diff(:,2));


%% TRAJECTORIES DATA FILE


nf = 2; % 1 

path2 = 'C:\Kezhi\MyCode!!!\ManualVideos\Check_Align_samples\';

% please add the folder name here
addpath(genpath([path2,'.']));

ffmpeg = 'C:\FFMPEG\bin\ffmpeg';
ffprobe = 'C:\FFMPEG\bin\ffprobe';


root = 'C:\Kezhi\MyCode!!!\ManualVideos\Check_Align_samples\';
hdf_folder = root;
trajectories_folder = [root,'Results\'];

root_folder = genpath([root,'.']);

traj_file=dir([trajectories_folder,'*_trajectories.hdf5']);
num_file = size(traj_file,1);

% for nf = 1:1;
% choose index of file

name  = traj_file(nf).name(1:end-18);
excel_name = [name, '.log.csv'];
hdf5_name = [name, '.hdf5'];

trajectories_file = [trajectories_folder,traj_file(nf).name];
plate_worms = h5read(trajectories_file, '/plate_worms');
timestamp = h5read(trajectories_file, '/timestamp/raw');
timestamp_time = h5read(trajectories_file, '/timestamp/time');

% import csv data
ss = importdata([root,excel_name]);

% read data: real time, media time, and stage coordinates xy
real_time = ss.textdata(2:end,1);
media_time = ss.textdata(2:end,2);
stage_xy = ss.data;

% check if number of rows in csv is equal to number of frames 
data_rows = size(ss.data,1);
if data_rows+1 ~= size(ss.textdata,1);
    error('excel file problem: number of rows are not in uniform')
end

% convert time from text to number 
media_time_vec = [];

% convert time to seconds value
for ii = 1:data_rows;
    str1 = media_time(ii);
    t1 = datevec(str1);
    media_time_vec(ii) = t1(5)*60 + t1(6);
end

hdf5_path = [hdf_folder,hdf5_name];


real_frame = timestamp(plate_worms.frame_number+1); %will match the indexes of segworm. Add one because the python indexing.
real_time_frame = timestamp_time(plate_worms.frame_number+1);


%% read key data
% mask information
mask_info = h5info(hdf5_path, '/mask');
% mask matrix, 3 dimentions
mask = h5read(hdf5_path, '/mask');
% size of each frame
frame_size = mask_info.Dataspace.Size(1:2);
frame_total = mask_info.Dataspace.Size(3);
frame_pos = h5read(hdf5_path,'/vid_frame_pos'); % start from 0
time_pos = h5read(hdf5_path,'/vid_time_pos');
normalize_val = 1000;

if frame_total ~= length(real_time_frame)
    error('number of frames is equal to number of time stamp in tragactories');
end

mask_central = [plate_worms.coord_x,plate_worms.coord_y,real_time_frame];

%% calculate the central moving difference between each frame, and save them in 'diff_mask_central(:,1,2)'
mask_bw = []; 
CC_max = logical(zeros(frame_size));
mask_previous = mask(:,:,1);

for ii = 1:frame_total-1;
    ii
    mask_backadjust = mask(:,:,ii);
    substitute_intensity = round(median(median(mask_backadjust(mask_backadjust>0)))*1.1);
    mask_backadjust(mask_backadjust<1)=substitute_intensity;

    mask_backadjust2 = mask(:,:,ii+1);
    mask_backadjust2(mask_backadjust2<1)=substitute_intensity;
    
%     level = graythresh(mask_backadjust);
%     mask_bw = im2bw(mask_backadjust,level)<0.5;
%     CC = bwconncomp(mask_bw);
    
%     numPixels = cellfun(@numel,CC.PixelIdxList);
%     [biggest,idx] = max(numPixels);
%     CC_max_this = CC_max;
%     CC_max_this(CC.PixelIdxList{1,idx}) = logical(1);
%     S = regionprops(CC_max_this,'Centroid');
%     mask_central(ii,1:2)=S.Centroid;

%     if ii == frame_total
%         break;
%     else
       pos_diff_pixel = double(mask_backadjust) - double(mask_backadjust2);
        sum_postive = sum(sum(pos_diff_pixel(pos_diff_pixel>20)));
        sum_negtive = sum(sum(pos_diff_pixel(pos_diff_pixel<-20)));
       pixel_diff(ii,1:2) = [sum_postive, sum_negtive]; 
%    end
end
% normalize time
%mask_central(:,3) = time_pos(1:size(mask_central,1))/normalize_val;
pixel_diff(:,3) = abs(pixel_diff(:,1).*pixel_diff(:,2));

%%

diff_mask_central = mask_central(2:end,1:2) - mask_central(1:end-1,1:2);
% delete outliers
diff_mask_central(abs(diff_mask_central(:,1))>30,1) = 0;
diff_mask_central(abs(diff_mask_central(:,2))>30,2) = 0;

diff_mask_central(:,3) = sqrt(diff_mask_central(:,1).^2+diff_mask_central(:,2).^2);
%diff_mask_central(:,4) = mask_central(2:end,3)/mask_central(end,3)*14.999*60;
diff_mask_central(:,4) = mask_central(2:end,3);

%% read stage moving information and save it in 'diff_mask_central(:,5,6)' according to time stamps in 'diff_mask_central(:,4)'
pt = 2;
for ii = 1:frame_total-1;
    ii
    if (pt <= data_rows) & (media_time_vec(pt)<diff_mask_central(ii,4))
          % switch x,y axis !!! (sometimes it is not necessary. Look out before use it)
        diff_mask_central(ii,5:1:6) = stage_xy(pt,:)- stage_xy(pt-1,:);
        pt = pt + 1;
    end 
end
num_pt = pt;

% calculate the absolute moving distance based on x-axis and y-axis
diff_mask_central(:,7) = sqrt(diff_mask_central(:,5).^2+diff_mask_central(:,6).^2);


% %% show some align figures
% figure, plot(diff_mask_central(1:10000,3)), hold on, plot(diff_mask_central(1:10000,7)/25,'red')
% figure, plot(diff_mask_central(10000:20000,3)), hold on, plot(diff_mask_central(10000:20000,7)/25,'red')
% figure, plot(diff_mask_central(end-10000:end,3)), hold on, plot(diff_mask_central(end-10000:end,7)/25,'red')

save_file = ['check_ex',num2str(nf),'.mat'];
save(save_file);

load_file = ['check_ex',num2str(nf),'.mat'];
load(load_file)

%% read pixels per microns from xml files
for qq = 1:2;
    if qq == 1
        findLabel1 = 'microns';
    elseif qq == 2
        findLabel1 = 'pixels';
    end

    % end - 8 because: excel name ends with '.log.csv', 8 letters
    xDoc1 = xmlread([root,excel_name(1:end-8),'.info.xml']);
    allListitems1 = xDoc1.getElementsByTagName(findLabel1);

    thisListitem1 = allListitems1.item(0);

    thisList1 = thisListitem1.getElementsByTagName('x');
    thisElement1 = thisList1.item(0);
    thisList2 = thisListitem1.getElementsByTagName('y');
    thisElement2 = thisList2.item(0);
    if qq == 1
        x_microns = str2num(thisElement1.getFirstChild.getData);
        y_microns = str2num(thisElement2.getFirstChild.getData);
    elseif qq ==2
        x_pixels = str2num(thisElement1.getFirstChild.getData);
        y_pixels = str2num(thisElement2.getFirstChild.getData);
    end
end
x_pixel_per_microns = x_pixels/x_microns;
y_pixel_per_microns = y_pixels/y_microns;

%% normalize

% normalize standard stage moving
diff_mask_central(:,5) = diff_mask_central(:,5)/y_pixel_per_microns;
diff_mask_central(:,6) = diff_mask_central(:,6)/x_pixel_per_microns;
diff_mask_central(:,7) = sqrt(diff_mask_central(:,5).^2+diff_mask_central(:,6).^2);
diff_mask_central(:,8) = pixel_diff(:,3)/10e8;

sorted_diff_area = sort(diff_mask_central(:,8),'descend');
% maximum keep ?% time slots that are potentially belongs to stage moving
pencentage_impulse = 0.05+(num_pt*6)/frame_total;
thresh_diff_area = sorted_diff_area(round(pencentage_impulse*size(sorted_diff_area,1)));
% fix 0 in diff_mask_central(:,8)
diff_mask_central((diff_mask_central(:,8)==0),8)=thresh_diff_area*2;
% threshold of difference of mask move
thresh_diff_mask = 3;

%% build a match
diff_leng = size(diff_mask_central,1);
moving_frame = zeros(diff_leng,2);
% identify moving frames: central moves >3 && areas change > threshold
moving_frame(:,1) = (diff_mask_central(:,3)>thresh_diff_mask)&(diff_mask_central(:,8)>thresh_diff_area);
moving_frame(:,2) = moving_frame(:,1);

% reduce each peak to "single frame length"
for qq = 1: diff_leng-2;
    if moving_frame(qq,2) == 1
        if moving_frame(qq+2,2) ==1
            moving_frame(max(1,qq-3):min(qq+1,diff_leng),2)=0;
            moving_frame(max(1,qq+3):min(qq+7,diff_leng),2)=0;
        else 
            moving_frame(max(1,qq-5):min(qq-1,diff_leng),2)=0;
            moving_frame(max(1,qq+1):min(qq+5,diff_leng),2)=0;
        end
    end
end

mov_fra_ind = (find(moving_frame(:,2)>0));
csv_ind = (find(diff_mask_central(:,7)>0));

% adjust the scale of "mov_fra_ind" to make it match to the "csv_ind"
match_mtx = zeros(20,40);
min_match_mtx_elem = 1e10;
min_fra_ind_match = [];
for mm1 = 1:20;
    for mm2 = 1:80;
        [match_mtx(mm1, mm2),mov_fra_ind_match] = cal_match_score(csv_ind, mov_fra_ind, mm1, 0.8+0.005*mm2, diff_mask_central(:,7), diff_mask_central(:,3));
        if match_mtx(mm1, mm2) < min_match_mtx_elem
            min_match_mtx_elem = match_mtx(mm1, mm2);
            min_fra_ind_match = mov_fra_ind_match;
            mm1_min_ind = mm1;
            mm2_min_ind = mm2;
        end
    end
end


%% refine aligment indexes

% gaussian window covers the impluse peak
% gaussian window size
filter_length = 81;
% half of the window size
half_filter = (filter_length-1)/2 ;
% generate the gausssian window vector
guass_window = gausswin(filter_length);

% parameters to calculate the 'shift_to_left'
shift_weights = exp(-2.4:0.6:0)';
normal_shift_weights = shift_weights/sum(shift_weights);
% only conder nearest 10 peaks to calculate the 'shift_to_left'
shift_consider = zeros(1,length(shift_weights));

mov_fra_ind_extend = round((mov_fra_ind( mm1_min_ind:end)-mov_fra_ind(mm1_min_ind))*(0.8+0.005*mm2_min_ind) + csv_ind(1));
mov_compare = zeros(diff_leng,1);
len_mov_fra_ind_extend = length(mov_fra_ind_extend);
for ll = 1:len_mov_fra_ind_extend
    mov_compare(mov_fra_ind_extend(ll)) = sum(diff_mask_central(max(1,mov_fra_ind(ll)-4):min(diff_leng,mov_fra_ind(ll)+4),3));
end

% absolute difference of stage central with extensions on both sides with
% all 0. 
abs_diff0 = abs([zeros(half_filter,1);mov_compare;zeros(half_filter,1)]);

% number of nonzeros peaks
no_nonzero_csv_ind = length(csv_ind);


% compensate the difference: shift impulse index to left 
shift_to_left = zeros(no_nonzero_csv_ind+1,1);

% determine x mask index
for ii = 1:no_nonzero_csv_ind;
    % 
    if ii == 425
        ii
    end
   %ind_consider0 = round((min_fra_ind_match(ii,1)-shift_to_left(ii)+min_fra_ind_match(ii,2))/2);
   ind_consider0 = min_fra_ind_match(ii,1)-shift_to_left(ii);
   
    % actually it should be: start_in_abs_diff = ind_consider -
    % half_filter+ half_filter; to find the indexes of values of interests
    % in 'abs_diff'
    start_in_abs_diff0 = max(1,ind_consider0-half_filter);
    % consider the case when 'end_inabs_diff0' is larger than 'diff_leng'
    end_in_abs_diff0 = min(diff_leng,ind_consider0+half_filter);
    % find the peak index with maximum weights, within the areas of 'positive/negtive areas are both large'
    sum_abs_diff0_sec = sum(start_in_abs_diff0+half_filter:end_in_abs_diff0+half_filter);
    peak_weights0 = min_fra_ind_match(ii,3)-abs_diff0(start_in_abs_diff0+half_filter:end_in_abs_diff0+half_filter).*guass_window(1:(end_in_abs_diff0-start_in_abs_diff0+1));
    % debug
    if sum_abs_diff0_sec==0
        start_in_abs_diff0
        end_in_abs_diff0
    end
    [weigths_max0, weights_max_ind0] = min(peak_weights0);
    % save results of indexes and absolute values of diff
    mask_ind(ii,1) = start_in_abs_diff0  + weights_max_ind0 - 1 ;
    mask_ind(ii,2) = sum(diff_mask_central(mask_ind(ii,1)-3:mask_ind(ii,1)+3,3));    
    
    % cancel this peak in abs_diff0
    abs_diff0(mask_ind(ii,1)+half_filter) = 1;
    
    % update 'shift_to_left'. it is caluclated based on last 10 shift
    % values
    shift_to_left(ii) = shift_to_left(ii) -weights_max_ind0 + half_filter;
    shift_consider(1:end-1) = shift_consider(2:end);
    %shift_consider(end) = min_fra_ind_match(ii,5)-mask_ind_x(ii,1);
    shift_consider(end) =  shift_to_left(ii);
    shift_to_left(ii+1) = round(shift_consider*normal_shift_weights);
end
align_extend_old = min_fra_ind_match(:,2);
min_fra_ind_match(:,2) = mask_ind(:,1);

%% important: obtain the alighment indexes here 
min_fra_ind_match(:,5) = round((min_fra_ind_match(:,2)-csv_ind(1))/(0.8+0.005*mm2_min_ind)+mov_fra_ind(mm1_min_ind));

% check result of this section
mov_fra_ind_extend = round((mov_fra_ind( mm1_min_ind:end)-mov_fra_ind(mm1_min_ind))*(0.8+0.005*mm2_min_ind) + csv_ind(1));
mov_fra_ind_Noextend = mov_fra_ind( mm1_min_ind:end) ;

mov_compare2 = zeros(diff_leng,1);

mov_compare2(min_fra_ind_match(:,5)) = 40;
figure, plot(diff_mask_central(:,7),'r-o');
hold on , 
plot(mov_compare,'y-');
plot(mov_compare2,'b-');
plot (diff_mask_central(:,3),'g-');

%% find corresponding impulses, for x
% initialize important values


% absolute difference of stage central with extensions on both sides with
% all 0. 
abs_diff = abs([zeros(half_filter,1);diff_mask_central(:,1);zeros(half_filter,1)]);
% possible absolute difference based on "difference multiplication"
possible_abs_diff =  [zeros(half_filter,1);(diff_mask_central(:,8)>thresh_diff_area);zeros(half_filter,1)];
% the indexes of nonzero peaks in csv
csv_ind_x =find( abs(diff_mask_central(:,5))>0);
% number of nonzeros peaks
no_nonzero_x = size(csv_ind_x,1);

diff_mask_central(:,9:10) = zeros(size(diff_mask_central,1),2);

% compensate the difference: shift impulse index to left 
shift_to_left = zeros(no_nonzero_x+1,1);
shift_consider = zeros(1,length(shift_weights));

% determine x mask index
for ii = 1:no_nonzero_x;
    % 
    %debug
    if ii == 35
        ii
    end
    ind_consider = min_fra_ind_match(ii,5)-shift_to_left(ii);
    % actually it should be: start_in_abs_diff = ind_consider -
    % half_filter+ half_filter; to find the indexes of values of interests
    % in 'abs_diff'
    start_in_abs_diff = ind_consider;
    end_in_abs_diff = start_in_abs_diff+filter_length-1;
    % find the peak index with maximum weights, within the areas of 'positive/negtive areas are both large'
    peak_weights = abs_diff(start_in_abs_diff:end_in_abs_diff).*guass_window.*possible_abs_diff(start_in_abs_diff:end_in_abs_diff);
    [weigths_max, weights_max_ind] = max(peak_weights);
    % save results of indexes and absolute values of diff
    mask_ind_x(ii,1) = start_in_abs_diff  + weights_max_ind - 1 -half_filter;
    mask_ind_x(ii,2) = diff_mask_central(mask_ind_x(ii,1),3);    
    
    % update 'shift_to_left'. it is caluclated based on last 10 shift
    % values
    shift_to_left(ii) = shift_to_left(ii) -weights_max_ind + half_filter;
    shift_consider(1:end-1) = shift_consider(2:end);
    %shift_consider(end) = min_fra_ind_match(ii,5)-mask_ind_x(ii,1);
    shift_consider(end) =  shift_to_left(ii);
    shift_to_left(ii+1) = round(shift_consider*normal_shift_weights);
end


%% find corresponding impulses, for y



% absolute difference of stage central with extensions on both sides with
% all 0. 
abs_diff = abs([zeros(half_filter,1);diff_mask_central(:,2);zeros(half_filter,1)]);
% the indexes of nonzero peaks in csv
csv_ind_y = find( abs(diff_mask_central(:,6))>0);
% number of nonzeros peaks
no_nonzero_y = size(csv_ind_y,1);

% compensate the difference: shift impulse index to left 
shift_to_left = zeros(no_nonzero_y+1,1);
shift_consider = zeros(1,length(shift_weights));

% determine y mask index
for ii = 1:no_nonzero_y;
    % 
    ind_consider = min_fra_ind_match(ii,5)-shift_to_left(ii);
    % actually it should be: start_in_abs_diff = ind_consider -
    % half_filter+ half_filter; to find the indexes of values of interests
    % in 'abs_diff'
    start_in_abs_diff = ind_consider;
    end_in_abs_diff = start_in_abs_diff+filter_length-1;
    % find the peak index with maximum weights, within the areas of 'positive/negtive areas are both large'
    peak_weights = abs_diff(start_in_abs_diff:end_in_abs_diff).*guass_window.*possible_abs_diff(start_in_abs_diff:end_in_abs_diff);
    [weigths_max, weights_max_ind] = max(peak_weights);
    % save results of indexes and absolute values of diff
    mask_ind_y(ii,1) = start_in_abs_diff  + weights_max_ind - 1 -half_filter;
    mask_ind_y(ii,2) = diff_mask_central(mask_ind_y(ii,1),3);    
    
    % update 'shift_to_left'. it is caluclated based on last 10 shift
    % values
    shift_consider(1:end-1) = shift_consider(2:end);
    shift_consider(end) = min_fra_ind_match(ii,5)-mask_ind_y(ii,1);
    shift_to_left(ii+1) = round(shift_consider*normal_shift_weights);
end


% figure, plot(diff_mask_central(1:1000,3)), hold on, plot(diff_mask_central(1:1000,7)/25,'red')
% %hold on, plot(mask_ind(1:1000,1),mask_ind(1:1000,2),'m*');
% hold on, plot(pixel_diff(1:1000,3)/10e9,'go-')


%% calculate stage moving vectors
% range1 is the left and right range of searching area, in alignment. So
% the the window size will be 2*range1+1
range1 = 15;

if mask_ind_x(1,1)< range1
     mask_ind_x = mask_ind_x(2:end,:);
end
if mask_ind_y(1,1)< range1
     mask_ind_y = mask_ind_y(2:end,:);
end
% initilaize results vector
stage_move_x = zeros(diff_leng,1);
stage_move_y = zeros(diff_leng,1);
% align csv_ind to mask_ind
for ii = 1:no_nonzero_x;
 %   ii
    % calculate stage moving vectors
    stage_move_x = align_func(stage_move_x, diff_mask_central(:,1), diff_mask_central(:,5),range1, mask_ind_x, csv_ind_x, ii);
end
for ii = 1:no_nonzero_y;
 %   ii
    % calculate stage moving vectors
    stage_move_y = align_func(stage_move_y, diff_mask_central(:,2), diff_mask_central(:,6),range1, mask_ind_y, csv_ind_y, ii);
end

diff_mask_central(:,9) = stage_move_x;
diff_mask_central(:,10) = stage_move_y;

% save result
temp_text_file = ['temp_text',num2str(nf),'.mat'];
save(temp_text_file, 'stage_move_x', 'stage_move_y','hdf5_path');

% end


%% test and draw
% 
% load([root,name,'_features.mat']);
% size(worm.posture.skeleton.x)
% figure,plot(worm.posture.skeleton.x(25, :), worm.posture.skeleton.y(25, :)); axis equal
% 
% res_x_central_diff = diff_mask_central(:,1)-diff_mask_central(:,9);
% res_y_central_diff = diff_mask_central(:,2)-diff_mask_central(:,10);
% res_x_central_diff(abs(res_x_central_diff)>4)=0;
% res_y_central_diff(abs(res_y_central_diff)>4)=0;
% 
% res_x_central = zeros(size(res_x_central_diff,1)+1,1); 
% res_y_central = zeros(size(res_y_central_diff,1)+1,1); 
% res_x_central(1) = worm.posture.skeleton.x(25, 1);
% res_y_central(1) = worm.posture.skeleton.y(25, 1);
% for ii = 2:size(res_x_central);
%     res_x_central(ii) = res_x_central(ii-1)+res_x_central_diff(ii-1);
%     res_y_central(ii) = res_y_central(ii-1)+res_y_central_diff(ii-1);
% end
% figure, plot(res_x_central, res_y_central); axis equal
% 
% % figure, plot(res_x_central(1), res_y_central(1),'r*'); axis equal
% % for ii= 2:length(res_x_central);
% %     hold on, plot(res_x_central(ii), res_y_central(ii),'r*');
% %     pause(1);
% % end
% 

