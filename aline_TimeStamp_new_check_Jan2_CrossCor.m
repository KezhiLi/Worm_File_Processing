clear
clc

%% TRAJECTORIES DATA FILE

% the index of file in consideration
nf = 1; % 1 

% the hdf5/csv/xml file folder
% path2 = 'C:\Kezhi\MyCode!!!\ManualVideos\Check_Align_samples\';
path2 = 'C:\Kezhi\MyCode!!!\ManualVideos\';

% make all subfolder available
addpath(genpath([path2,'.']));

ffmpeg = 'C:\FFMPEG\bin\ffmpeg';
ffprobe = 'C:\FFMPEG\bin\ffprobe';

% the root 
root = 'Z:\thecus\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\15-03-11\';
hdf_folder = 'Z:\MaskedVideos\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\15-03-11\';
trajectories_folder = 'Z:\Results\nas207-1\experimentBackup\from pc207-7\!worm_videos\copied_from_pc207-8\Andre\15-03-11\';

root_folder = genpath([root,'.']);

traj_file=dir([trajectories_folder,'*_trajectories.hdf5']);
num_file = size(traj_file,1);
good_ind = [];

save path2.mat
for nf = 1: length(traj_file);
%for nf = 1:length(traj_file);
% choose index of file

 try

name  = traj_file(nf).name(1:end-18);
excel_name = [name, '.log.csv'];
hdf5_name = [name, '.hdf5'];

% read info from tracjectories
trajectories_file = [trajectories_folder,traj_file(nf).name];
plate_worms = h5read(trajectories_file, '/plate_worms');
timestamp = h5read(trajectories_file, '/timestamp/raw');
timestamp_time = h5read(trajectories_file, '/timestamp/time');

% import csv data
%ss = importdata([root,excel_name]);
ss = read_csv_data(root,excel_name);

% read csv data: real time, media time, and stage coordinates xy
real_time = ss.textdata(2:end,1);
media_time = ss.textdata(2:end,2);
stage_xy = ss.data;

% check if number of rows in csv is equal to number of frames 
data_rows = size(ss.data,1);
if data_rows+1 ~= size(ss.textdata,1);
    error('excel file problem: number of rows are not in uniform')
end



% convert time to seconds value
if isa(media_time,'double')
    media_time_vec = media_time;
else
    % convert time from text to number 
    media_time_vec = size(media_time);
    for ii = 1:data_rows;
        str1 = media_time(ii);
        t1 = datevec(str1);
        media_time_vec(ii) = t1(5)*60 + t1(6);
    end
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
% calculate the difference of the central of the mask
diff_mask_central = mask_central(2:end,1:2) - mask_central(1:end-1,1:2);
% delete outliers
diff_mask_central(abs(diff_mask_central(:,1))>30,1) = 0;
diff_mask_central(abs(diff_mask_central(:,2))>30,2) = 0;

% column 3 is the absolute shift distance considering both x,y direcition
diff_mask_central(:,3) = sqrt(diff_mask_central(:,1).^2+diff_mask_central(:,2).^2);
%diff_mask_central(:,4) = mask_central(2:end,3)/mask_central(end,3)*14.999*60;
diff_mask_central(:,4) = mask_central(2:end,3);

%% read stage moving information and save it in 'diff_mask_central(:,5,6)' according to time stamps in 'diff_mask_central(:,4)'
pt = 2;
for ii = 1:frame_total-1;
    ii
    % insert the stage motion in the right time
    if (pt <= data_rows) & (media_time_vec(pt)<diff_mask_central(ii,4))
          % estimate the stage motion in cvs
        diff_mask_central(ii,5:1:6) = stage_xy(pt,:)- stage_xy(pt-1,:);
        % find next shift
        pt = pt + 1;
    end 
end
num_pt = pt;

% %% show some align figures
% figure, plot(diff_mask_central(1:10000,3)), hold on, plot(diff_mask_central(1:10000,7)/25,'red')
% figure, plot(diff_mask_central(10000:20000,3)), hold on, plot(diff_mask_central(10000:20000,7)/25,'red')
% figure, plot(diff_mask_central(end-10000:end,3)), hold on, plot(diff_mask_central(end-10000:end,7)/25,'red')

% save_file = ['check_ex',num2str(nf),'.mat'];
% save(save_file);
% 
% load_file = ['check_ex',num2str(nf),'.mat'];
% load(load_file)

% save_file = [name,'-data.mat'];
% save(save_file);
% 
% load_file = [name,'-data.mat'];
% load(load_file)

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
moving_frame(:,1) = (((abs(xShift)>4)|(abs(yShift)>4))>0)&(diff_mask_central(:,3)>1.5);
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
    for mm2 = 1:40;
        [match_mtx(mm1, mm2),mov_fra_ind_match] = cal_match_score(csv_ind, mov_fra_ind, mm1, 0.9+0.005*mm2, diff_mask_central(:,7), diff_mask_central(:,3));
        if match_mtx(mm1, mm2) < min_match_mtx_elem
            min_match_mtx_elem = match_mtx(mm1, mm2);
            min_fra_ind_match = mov_fra_ind_match;
            mm1_min_ind = mm1;
            mm2_min_ind = mm2;
        end
    end
end
% calculate the orignal indexes
min_fra_ind_match(:,5) = round((min_fra_ind_match(:,2)-csv_ind(1))/(0.9+0.005*mm2_min_ind)+mov_fra_ind(mm1_min_ind));

%% refine aligment indexes

% gaussian window covers the impluse peak
% gaussian window size
filter_length = 71;
% half of the window size
half_filter = (filter_length-1)/2 ;
% generate the gausssian window vector
guass_window = gausswin(filter_length);

% parameters to calculate the 'shift_to_left'
shift_weights = exp(-2.4:0.6:0)';
normal_shift_weights = shift_weights/sum(shift_weights);
% only conder nearest 10 peaks to calculate the 'shift_to_left'
shift_consider = zeros(1,length(shift_weights));

mov_compare = zeros(diff_leng,2);
len_mov_fra_ind_extend = length(mov_fra_ind);
for ll = 1:len_mov_fra_ind_extend
    mov_compare(mov_fra_ind(ll),1:2) = sum(diff_mask_central(max(1,mov_fra_ind(ll)-3):min(diff_leng,mov_fra_ind(ll)+3),1:2));
end

% absolute difference of stage central with extensions on both sides with
% all 0. 
diff0 = [zeros(half_filter,2);mov_compare;zeros(half_filter,2)];

% number of nonzeros peaks
no_nonzero_csv_ind = length(csv_ind);


% compensate the difference: shift impulse index to left 
shift_to_left = zeros(no_nonzero_csv_ind+1,1);

mask_ind = zeros(no_nonzero_csv_ind,2);
% determine x mask index
for ii = 1:no_nonzero_csv_ind;
    % debug purpose
%     if ii == 93
%         ii
%     end
  
   ind_consider0 = round((min_fra_ind_match(ii,1)+min_fra_ind_match(ii,5))/2)-shift_to_left(ii);
   
    % actually it should be: start_in_abs_diff = ind_consider -
    % half_filter+ half_filter; to find the indexes of values of interests
    % in 'abs_diff'
    start_in_abs_diff0 = max(1,ind_consider0-half_filter);
    % consider the case when 'end_inabs_diff0' is larger than 'diff_leng'
    end_in_abs_diff0 = min(diff_leng,ind_consider0+half_filter);
    % find the peak index with maximum weights, within the areas of 'positive/negtive areas are both large'
    sum_abs_diff0_sec = sum(ind_consider0:end_in_abs_diff0+half_filter);
    
    ind_in_consider = ind_consider0:end_in_abs_diff0+half_filter;
    % use a guassian window to esitmate the best weights
    peak_weights0 = ones(length(ind_in_consider),1)*diff_mask_central(min_fra_ind_match(ii,1),5:6)-diff0(ind_in_consider,1:2).*(guass_window(round((filter_length-length(ind_in_consider))/2)+(1:length(ind_in_consider)))*ones(1,2));
    % debug use
    if sum_abs_diff0_sec==0
        start_in_abs_diff0
        end_in_abs_diff0
        error('no impluses in the interval')
    end
    [weigths_max0, weights_max_ind0] = min(sum(abs(peak_weights0),2));
    % save results of indexes and absolute values of diff
    mask_ind(ii,1) = ind_in_consider(weights_max_ind0)-half_filter ;
    mask_ind(ii,2) = sum(diff_mask_central(max(1,mask_ind(ii,1)-3):min(diff_leng,mask_ind(ii,1)+3),3));    
    
    % cancel this peak in abs_diff0
    diff0(mask_ind(ii,1)+half_filter,1:2) = 1;
    
    % update 'shift_to_left'. it is caluclated based on last 10 shift
    % values
    shift_to_left(ii) = shift_to_left(ii) -weights_max_ind0 + half_filter;
    shift_consider(1:end-1) = shift_consider(2:end);
    %shift_consider(end) = min_fra_ind_match(ii,5)-mask_ind_x(ii,1);
    shift_consider(end) =  shift_to_left(ii);
    shift_to_left(ii+1) = round(shift_consider*normal_shift_weights);
end
align_extend_old = min_fra_ind_match(:,2);
% use the sorted indexes as results
min_fra_ind_match(:,6) = sort(mask_ind(:,1));

%% choose better result from 'min_fra_ind_match(:,5)' and 'min_fra_ind_match(:,9)'
shift_para = 1;
min_fra_ind_match(:,7) = zeros(size(min_fra_ind_match(:,6)));
for iin = 1: no_nonzero_csv_ind;
    if iin == 93
        iin
    end
    % compare two alignment results, if they are the same, then choose any
    % one of them
   if (min_fra_ind_match(iin,5) ==  min_fra_ind_match(iin,6)) | (iin == 1)
       min_fra_ind_match(iin,7) = min_fra_ind_match(iin,5);
   else
       % if two alignments are not the same, compare the shift in x,y axis,
       % in terms of 'shift of central of area' and 'result of cross
       % correlation', respectively. 
       
       % magnitudes of shift of area centre
       xy_sum_diffCentr_5 = sum(diff_mask_central(max(1,min_fra_ind_match(iin,5)-3):...
           min(diff_leng,min_fra_ind_match(iin,5)+3),1:2),1);
       % magnitudes of shift of area centre
       xy_sum_diffCentr_6 = sum(diff_mask_central(max(1,min_fra_ind_match(iin,6)-3):...
           min(diff_leng,min_fra_ind_match(iin,6)+3),1:2),1);
       % magnitudes of the sum of crossCorrelation around each peak, for
       % the 5th column of 'min_fra_ind_match'
       xy_sum_CrossCor_5 = [sum(xShift(max(1,min_fra_ind_match(iin,5)-3):...
           min(diff_leng,min_fra_ind_match(iin,5)+3))),...
           sum(yShift(max(1,min_fra_ind_match(iin,5)-3):...
           min(diff_leng,min_fra_ind_match(iin,5)+3)))];
       % magnitudes of the sum of crossCorrelation around each peak, for the 6th column of 'min_fra_ind_match' 
       xy_sum_CrossCor_6 = [sum(xShift(max(1,min_fra_ind_match(iin,6)-3):...
           min(diff_leng,min_fra_ind_match(iin,6)+3))),...
           sum(yShift(max(1,min_fra_ind_match(iin,6)-3):...
           min(diff_leng,min_fra_ind_match(iin,6)+3)))];
       xy_csv = diff_mask_central(min_fra_ind_match(iin,1),5:6);
       % estimate the error in terms of col
       error_5 = sum(abs(xy_sum_diffCentr_5 - xy_csv)+abs(xy_sum_CrossCor_5 - xy_csv))...
           + shift_para*abs((min_fra_ind_match(iin,1)-min_fra_ind_match(iin-1,1))...
             -(min_fra_ind_match(iin,5)-min_fra_ind_match(iin-1,7)));
       error_6 = sum(abs(xy_sum_diffCentr_6 - xy_csv)+abs(xy_sum_CrossCor_6 - xy_csv))...
           + shift_para*abs((min_fra_ind_match(iin,1)-min_fra_ind_match(iin-1,1))...
             -(min_fra_ind_match(iin,6)-min_fra_ind_match(iin-1,7)));
       
       % select the alignment result with smaller x,y shift errors and
       % non-repeated indexes. 
       if (error_5 <= error_6) & ~any(abs(min_fra_ind_match(1:iin,7)-min_fra_ind_match(iin,5))<1e-5)
            min_fra_ind_match(iin,7) = min_fra_ind_match(iin,5);
       elseif (error_6 < error_5) & ~any(abs(min_fra_ind_match(1:iin,7)-min_fra_ind_match(iin,6))<1e-5)
            min_fra_ind_match(iin,7) = min_fra_ind_match(iin,6);
       else
           back_ind = 1;
           min_fra_ind_match(iin,7) = min_fra_ind_match(iin,5);
           % if there is a duplicated index, choose the 5th column as the
           % result
           while length(unique(min_fra_ind_match(1:iin,7)))~=iin
                min_fra_ind_match(iin-back_ind,7)=min_fra_ind_match(iin-back_ind,5);
                % the backwards index, increase from 1, until no duplicated
                % index exist in final index result vector
                back_ind = back_ind+1;
           end
       end
   end
   
end



%% important: obtain the alighment indexes here 

% check result of this section
mov_fra_ind_extend = round((mov_fra_ind( mm1_min_ind:end)-mov_fra_ind(mm1_min_ind))*(0.8+0.005*mm2_min_ind) + csv_ind(1));
mov_fra_ind_Noextend = mov_fra_ind( mm1_min_ind:end) ;

mov_compare2 = zeros(diff_leng,1);

mov_compare2(min_fra_ind_match(:,7)) = 35;
figure, plot(diff_mask_central(:,7),'r-o');
hold on , 
plot(sqrt(mov_compare(:,1).^2+mov_compare(:,1).^2),'y-');
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
temp_text_file = [name,'_align.mat'];
save([root,temp_text_file], 'nf','diff_mask_central', 'min_fra_ind_match','hdf5_path');

mov_frame_compare = [diff_mask_central(:,1), xShift, diff_mask_central(:,5),...
    diff_mask_central(:,2), yShift, diff_mask_central(:,6)];
% end

%% show all skeletons

% r
skeleton_hdf5 = h5read([trajectories_folder,name,'_skeletons.hdf5'],'/skeleton');

x_ske = (reshape(skeleton_hdf5(1,:,:), size(skeleton_hdf5,2),size(skeleton_hdf5,3)))';
y_ske = (reshape(skeleton_hdf5(2,:,:), size(skeleton_hdf5,2),size(skeleton_hdf5,3)))';

stage_mov_x_cum = cumsum(diff_mask_central(:,9));
stage_mov_y_cum = cumsum(diff_mask_central(:,10));

if size(x_ske,1) == length(stage_mov_x_cum)+1
    x_ske_cum = x_ske(2:end,:) + stage_mov_x_cum*ones(1,size(x_ske,2));
    y_ske_cum = y_ske(2:end,:) + stage_mov_y_cum*ones(1,size(y_ske,2));
else
    error('frame number size does not match');
end

show_fra_ind = round(size(x_ske_cum,1));
half_x = 320;
half_y = 240;

x_ske_cum_adjusted = x_ske_cum - min(x_ske_cum(1:show_fra_ind))+half_x*2;
y_ske_cum_adjusted = y_ske_cum - min(y_ske_cum(1:show_fra_ind))+half_y*2;

fig63 = figure(63),plot(x_ske_cum_adjusted(1,:),y_ske_cum_adjusted(1,:) ); axis equal; hold on
for tt_1 = 2:20:show_fra_ind;
    plot(x_ske_cum_adjusted(tt_1,:),y_ske_cum_adjusted(tt_1,:) );
end
hold off,
%% show and save as video

% initialize current image
curr_img = uint8(zeros(1800,1800));
% half size of window in focus
half_x = 320;
half_y = 240;

x_centr = -diff_mask_central(:,9);
y_centr = -diff_mask_central(:,10);
x_centr(1) = x_centr(1) + skeleton_hdf5(1,25, 1);
y_centr(1) = y_centr(1) + skeleton_hdf5(2,25, 1);

% x_centr_summ = round(cumsum(x_centr)*abs(x_pixel_per_microns));
% y_centr_summ = round(cumsum(y_centr)*abs(y_pixel_per_microns));
x_centr_summ = round(cumsum(x_centr));
y_centr_summ = round(cumsum(y_centr));

show_fra_ind = 1500;
x_centr_summ_adjusted = x_centr_summ - min(x_centr_summ(1:show_fra_ind))+half_x*2;
y_centr_summ_adjusted = y_centr_summ - min(y_centr_summ(1:show_fra_ind))+half_y*2;



for ii = 1:show_fra_ind;
    ii
    mask_backadjust = mask(:,:,ii+1);
    substitute_intensity = round(median(median(mask_backadjust(mask_backadjust>0)))*1.1);
    mask_backadjust(mask_backadjust<1)=substitute_intensity;
    
    mask_backadjust(1:4,:) = 256;
    mask_backadjust(end-3:end,:) = 256;
    mask_backadjust(:,1:4) = 256;
    mask_backadjust(:,end-3:end) = 256;
    
    curr_img = uint8(zeros(2048,2048));
    curr_img((y_centr_summ_adjusted(ii)-half_y):(y_centr_summ_adjusted(ii)+half_y-1),...
       (x_centr_summ_adjusted(ii)-half_x):(x_centr_summ_adjusted(ii)+half_x-1) ) = (mask_backadjust)'; 
    curr_img = curr_img(1:4:end, 1:4:end);
    curr_img(1:32:end,1:32:end)= 256;
%     curr_img(2:16:end,1:16:end)= 256;
%     curr_img(1:16:end,2:16:end)= 256;
%     curr_img(2:16:end,2:16:end)= 256;
    
    figure(10+nf), imshow(curr_img); 
    fps = 30;
      mov(ii) = save_crt_fra(name,ii, fps);
end

fname = [trajectories_folder,name,'.avi' ];
% movie2avi(mov, fname, 'compression', 'MSVC', 'fps', fps);

myVideo = VideoWriter(fname);
myVideo.Quality = 50;    % Default 75
open(myVideo);
writeVideo(myVideo, mov);
close(myVideo);

   good_ind = [good_ind, nf];
   save('path2_good_ind.mat', 'good_ind');
   %saveas(figure(1), 'testfig.fig');
   savefig(fig63,[trajectories_folder,name,'-ske.fig']);
   
  clear 
  load path2.mat
  load path2_good_ind.mat 
  
%   
catch
    continue;
end
    


% ske_traj_hdf5 = h5read([trajectories_folder,name,'_skeletons.hdf5'],'/trajectories_data');
% 
% ind_ii = 1;
% figure(51),plot(skeleton_hdf5(1,:,ind_ii), skeleton_hdf5(2,:,ind_ii)); axis equal; hold on
% for ind_ii = 101:100:size(skeleton_hdf5,3);
%     plot(skeleton_hdf5(1,:, ind_ii), skeleton_hdf5(2,:, ind_ii)); axis equal; hold on
% end
% 
% ind_ii = 1;
% figure(52),plot(plate_worms.coord_x(ind_ii), skeleton_hdf5(2,:,ind_ii)); axis equal; hold on
% for ind_ii = 101:100:size(skeleton_hdf5,3);
%     plot(skeleton_hdf5(1,:, ind_ii), skeleton_hdf5(2,:, ind_ii)); axis equal; hold on
% end
% 
% % xxx = skeleton_hdf5(1,25,1:10:end);
% % yyy = skeleton_hdf5(2,25,1:10:end);
% % figure(61),plot(reshape(xxx,size(xxx,3),1),reshape(yyy,size(yyy,3),1) ); axis equal;
% 
% cood_xx = ske_traj_hdf5.coord_x(1:10:end);
% cood_yy = ske_traj_hdf5.coord_y(1:10:end);
% figure(62),plot(cood_xx,cood_yy); axis equal;


end


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

