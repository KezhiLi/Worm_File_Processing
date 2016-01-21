clear
clc

%% TRAJECTORIES DATA FILE

% the index of file in consideration
nf = 3; % 1 

% the hdf5/csv/xml file folder
path2 = 'C:\Kezhi\MyCode!!!\ManualVideos\Check_Align_samples\';

% make all subfolder available
addpath(genpath([path2,'.']));

ffmpeg = 'C:\FFMPEG\bin\ffmpeg';
ffprobe = 'C:\FFMPEG\bin\ffprobe';

% the root 
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

% read info from tracjectories
trajectories_file = [trajectories_folder,traj_file(nf).name];
plate_worms = h5read(trajectories_file, '/plate_worms');
timestamp = h5read(trajectories_file, '/timestamp/raw');
timestamp_time = h5read(trajectories_file, '/timestamp/time');

% import csv data
ss = importdata([root,excel_name]);

% read csv data: real time, media time, and stage coordinates xy
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

% hdf5 file path
hdf5_path = [hdf_folder,hdf5_name];

% calculate the real time according to trajactories
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
% 
% for ii = 1:frame_total-1;
%     ii
%     mask_backadjust = mask(:,:,ii);
%     substitute_intensity = round(median(median(mask_backadjust(mask_backadjust>0)))*1.1);
%     mask_backadjust(mask_backadjust<1)=substitute_intensity;
% 
%     mask_backadjust2 = mask(:,:,ii+1);
%     mask_backadjust2(mask_backadjust2<1)=substitute_intensity;
%     
%     % the multiplication of difference of 2 frames
%        pos_diff_pixel = double(mask_backadjust) - double(mask_backadjust2);
%         sum_postive = sum(sum(pos_diff_pixel(pos_diff_pixel>20)));
%         sum_negtive = sum(sum(pos_diff_pixel(pos_diff_pixel<-20)));
%        pixel_diff(ii,1:2) = [sum_postive, sum_negtive]; 
% %    end
% end
% % normalize time
% %mask_central(:,3) = time_pos(1:size(mask_central,1))/normalize_val;
% pixel_diff(:,3) = abs(pixel_diff(:,1).*pixel_diff(:,2));

%% calculate shift from cross correlatoin
% set parameters
timeDiff = 1; % how many frames between aligned images?
dS = 2; % pixel downsampling factor (2 means half size)

% estimate transformation from one image frame to another
No_mask = size(mask, 3);
xShift = NaN(No_mask-timeDiff, 1);
yShift = NaN(No_mask-timeDiff, 1);
for ii = 1+timeDiff:No_mask
    % show the processing percentage
    disp(ii/No_mask)
    
    % subsample the image in before frame and after frame
    frame_bef = mask(1:dS:end, 1:dS:end, ii);
    frame_aft = mask(1:dS:end, 1:dS:end, ii - timeDiff);
    
    % use sum to find the worms in 2 frame, and the background of value 0
    frame_sum = abs(frame_bef)+abs(frame_aft)>0;
    frm_sum_col = sum(frame_sum,1);
    frm_sum_row = sum(frame_sum,2);
    
    % find the joint worm body
    worm_ind_col = find(frm_sum_col>0);
    worm_ind_row = find(frm_sum_row>0);
    
    % the pixel buffer round the worm body area
    pix_buffer = 5;
    
    % create a square of pixels that cover the worm body
    frame_bef = frame_bef(max(1,worm_ind_row(1)-pix_buffer)...
        :min(length(frm_sum_row),worm_ind_row(end)+pix_buffer),...
    max(1,worm_ind_col(1)-pix_buffer)...
        :min(length(frm_sum_col),worm_ind_col(end)+pix_buffer));
    frame_aft = frame_aft(max(1,worm_ind_row(1)-pix_buffer)...
        :min(length(frm_sum_row),worm_ind_row(end)+pix_buffer),...
    max(1,worm_ind_col(1)-pix_buffer)...
        :min(length(frm_sum_col),worm_ind_col(end)+pix_buffer));
    
    % estimate shift between images
    transMat = imregcorr(frame_bef , ...
       frame_aft  , 'translation');
   
    % calculate the shift in x,y directions
    xShift(ii - timeDiff) = transMat.T(3, 1)*dS;
    yShift(ii - timeDiff) = transMat.T(3, 2)*dS;
    
    xShift_temp = xShift;
    xShift = -yShift;
    yShift = -xShift_temp;
    
    
end


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
% please be care about the y_pixel and x_pixle here!!!
diff_mask_central(:,5) = diff_mask_central(:,5)/x_pixel_per_microns;
diff_mask_central(:,6) = diff_mask_central(:,6)/y_pixel_per_microns;
diff_mask_central(:,7) = sqrt(diff_mask_central(:,5).^2+diff_mask_central(:,6).^2);
% diff_mask_central(:,8) = pixel_diff(:,3)/10e8;
diff_mask_central(:,8) = diff_mask_central(:,3);

% 
% sorted_diff_area = sort(diff_mask_central(:,8),'descend');
% % maximum keep ?% time slots that are potentially belongs to stage moving
% pencentage_impulse = 0.05+(num_pt*6)/frame_total;
% thresh_diff_area = sorted_diff_area(round(pencentage_impulse*size(sorted_diff_area,1)));
% % fix 0 in diff_mask_central(:,8)
% diff_mask_central((diff_mask_central(:,8)==0),8)=thresh_diff_area*2;
% % threshold of difference of mask move
% thresh_diff_mask = 3;

%% build a match
diff_leng = size(diff_mask_central,1);
moving_frame = zeros(diff_leng,2);
% % identify moving frames: central moves >3 && areas change > threshold
% moving_frame(:,1) = (diff_mask_central(:,3)>thresh_diff_mask)&(diff_mask_central(:,8)>thresh_diff_area);
% moving_frame(:,2) = moving_frame(:,1);
moving_frame(:,1) = ((abs(xShift)>2)|(abs(yShift)>2))>0;
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
% calculate the orignal indexes
min_fra_ind_match(:,5) = round((min_fra_ind_match(:,2)-csv_ind(1))/(0.8+0.005*mm2_min_ind)+mov_fra_ind(mm1_min_ind));

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

% mov_fra_ind_extend = round((mov_fra_ind( mm1_min_ind:end)-mov_fra_ind(mm1_min_ind))*(0.8+0.005*mm2_min_ind) + csv_ind(1));
mov_compare = zeros(diff_leng,1);
len_mov_fra_ind_extend = length(mov_fra_ind);
for ll = 1:len_mov_fra_ind_extend
    mov_compare(mov_fra_ind(ll)) = sum(diff_mask_central(max(1,mov_fra_ind(ll)-4):min(diff_leng,mov_fra_ind(ll)+4),3));
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
    if ii == 49
        ii
    end
  
   ind_consider0 = round((min_fra_ind_match(ii,1)+min_fra_ind_match(ii,5))/2)-shift_to_left(ii);
   
    % actually it should be: start_in_abs_diff = ind_consider -
    % half_filter+ half_filter; to find the indexes of values of interests
    % in 'abs_diff'
    start_in_abs_diff0 = max(1,ind_consider0-half_filter);
    % consider the case when 'end_inabs_diff0' is larger than 'diff_leng'
    end_in_abs_diff0 = min(diff_leng,ind_consider0+half_filter);
    % find the peak index with maximum weights, within the areas of 'positive/negtive areas are both large'
    sum_abs_diff0_sec = sum(start_in_abs_diff0+half_filter:end_in_abs_diff0+half_filter);
    peak_weights0 = min_fra_ind_match(ii,3)-abs_diff0(start_in_abs_diff0+half_filter:end_in_abs_diff0+half_filter).*guass_window(1:(end_in_abs_diff0-start_in_abs_diff0+1));
    % debug use
    if sum_abs_diff0_sec==0
        start_in_abs_diff0
        end_in_abs_diff0
        error('no impluses in the interval')
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
min_fra_ind_match(:,6) = mask_ind(:,1);

%% choose better result from 'min_fra_ind_match(:,5)' and 'min_fra_ind_match(:,9)'

for iin = 1: no_nonzero_csv_ind;
    % compare two alignment results, if they are the same, then choose any
    % one of them
   if min_fra_ind_match(iin,5) ==  min_fra_ind_match(iin,6)
       min_fra_ind_match(iin,7) = min_fra_ind_match(iin,5);
   else
       % if two alignments are not the same, compare the shift in x,y axis,
       % in terms of 'shift of central of area' and 'result of cross
       % correlation', respectively. 
       xy_sum_diffCentr_5 = sum(diff_mask_central(max(1,min_fra_ind_match(iin,5)-3):...
           min(diff_leng,min_fra_ind_match(iin,5)+3),1:2),1);
       xy_sum_diffCentr_6 = sum(diff_mask_central(max(1,min_fra_ind_match(iin,6)-3):...
           min(diff_leng,min_fra_ind_match(iin,6)+3),1:2),1);
       xy_sum_CrossCor_5 = [sum(xShift(max(1,min_fra_ind_match(iin,5)-3):...
           min(diff_leng,min_fra_ind_match(iin,5)+3))),...
           sum(yShift(max(1,min_fra_ind_match(iin,5)-3):...
           min(diff_leng,min_fra_ind_match(iin,5)+3)))];
       xy_sum_CrossCor_6 = [sum(xShift(max(1,min_fra_ind_match(iin,6)-3):...
           min(diff_leng,min_fra_ind_match(iin,6)+3))),...
           sum(yShift(max(1,min_fra_ind_match(iin,6)-3):...
           min(diff_leng,min_fra_ind_match(iin,6)+3)))];
       xy_csv = diff_mask_central(min_fra_ind_match(iin,1),5:6);
       error_5 = sum(abs(xy_sum_diffCentr_5 - xy_csv)+abs(xy_sum_CrossCor_5 - xy_csv));
       error_6 = sum(abs(xy_sum_diffCentr_6 - xy_csv)+abs(xy_sum_CrossCor_6 - xy_csv));
       
       % select the alignment result with smaller x,y shift errors and
       % non-repeated indexes. 
       if (error_5 <= error_6) & ~any(abs(min_fra_ind_match(1:iin,7)-min_fra_ind_match(iin,5))<1e-5)
            min_fra_ind_match(iin,7) = min_fra_ind_match(iin,5);
       elseif (error_6 < error_5) & ~any(abs(min_fra_ind_match(1:iin,7)-min_fra_ind_match(iin,5))<1e-5)
            min_fra_ind_match(iin,7) = min_fra_ind_match(iin,6);
       else
           error('alignment error code 2');
       end
   end
end



%% important: obtain the alighment indexes here 

% check result of this section
mov_fra_ind_extend = round((mov_fra_ind( mm1_min_ind:end)-mov_fra_ind(mm1_min_ind))*(0.8+0.005*mm2_min_ind) + csv_ind(1));
mov_fra_ind_Noextend = mov_fra_ind( mm1_min_ind:end) ;

mov_compare2 = zeros(diff_leng,1);

mov_compare2(min_fra_ind_match(:,6)) = 35;
figure, plot(diff_mask_central(:,7),'r-o');
hold on , 
plot(mov_compare,'y-');
plot(mov_compare2,'b-');
plot (diff_mask_central(:,3),'g-');

stage_move_x = zeros(diff_leng,1);
stage_move_y = zeros(diff_leng,1);

range1 = 15;
% csv_ind_x =find( abs(diff_mask_central(:,5))>0);
% csv_ind_y =find( abs(diff_mask_central(:,6))>0);
% no_nonzero_x = size(csv_ind_x,1);
% no_nonzero_y = size(csv_ind_y,1);

% align csv_ind to mask_ind
for ii = 1:no_nonzero_csv_ind;
 %   ii
    % calculate stage moving vectors
    stage_move_x = align_func(stage_move_x, diff_mask_central(:,1), diff_mask_central(:,5),range1,  min_fra_ind_match(:,7), csv_ind, ii);
end
for ii = 1:no_nonzero_csv_ind;
 %   ii
    % calculate stage moving vectors
    stage_move_y = align_func(stage_move_y, diff_mask_central(:,2), diff_mask_central(:,6),range1,  min_fra_ind_match(:,7), csv_ind, ii);
end

diff_mask_central(:,9) = stage_move_x;
diff_mask_central(:,10) = stage_move_y;

% save result
temp_text_file = ['temp_text',num2str(nf),'.mat'];
save(temp_text_file, 'stage_move_x', 'stage_move_y','hdf5_path');

mov_frame_compare = [diff_mask_central(:,1), xShift, diff_mask_central(:,5),...
    diff_mask_central(:,2), yShift, diff_mask_central(:,6)];
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

