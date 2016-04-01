function [stage_vec, pixel_to_micro] = Align_TimeStamp_Func_ver2(masked_image_file, skeletons_file, features_mat, is_swimming)
%% read pixels per microns attribute. This info was extracted from the info.xml file
x_pixel_per_microns = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_x');
y_pixel_per_microns = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_y');

%% read stage data from the mask hdf5 file. This information was extracted
%from the .log.csv file 
stage_data = h5read(masked_image_file, '/stage_data');
%correct for duplicated data keeping the last instance of a given time
[stage_time, ind] = unique(stage_data.stage_time, 'last');
stage_time = stage_time*60; %Convert the time in seconds. It was originally in minutes.
stage_xy = [stage_data.stage_x(ind) , stage_data.stage_y(ind)];


%% the centroid position and real time frame it's better stored in the skeletons_file

video_timestamp_ind = h5read(skeletons_file, '/timestamp/raw');
video_timestamp_time = h5read(skeletons_file, '/timestamp/time');
trajectories_data = h5read(skeletons_file, '/trajectories_data');

%At this point after the skeletons analysis there should be only one valid
%index.
assert(all(trajectories_data.worm_index_joined==1))

real_time_frame = trajectories_data.timestamp_time;
mask_central = [trajectories_data.cnt_coord_x, trajectories_data.cnt_coord_y, real_time_frame];

%% test the timestamp make sense
% if the video is too short, report it and jump to next video
if real_time_frame(end)< 15
    error('the video is too short, less than 15 seconds: \n timestamp_time(end)=%d ',video_timestamp_time(end));
end


%% calculate shift from cross correlation between frames, and get the absolute difference between images
[xShift, yShift, abs_diff_fra] = shiftCrossCorrelation(masked_image_file);

%% get possible peaks from the pixel differences
% threshold(otsu) and trim 'abs_diff_fra' to make it good for tell the peaks
max_abs_diff = max(abs_diff_fra);
thresh = (graythresh(abs_diff_fra/max_abs_diff)*max_abs_diff)*0.8; % *0.8 to have a loose threshold

abs_diff_fra_th = abs_diff_fra;
abs_diff_fra_th(abs_diff_fra<thresh) = 0;   
abs_diff_fra_trim = imerode(abs_diff_fra_th, [1;1;1]);

%% main matrix used in alignment
diff_mask_central = zeros(size(mask_central,1)-1, 8);

% column 1  and 2 is the difference of the central of the mask in x and y
diff_mask_central(:, 1:2) = mask_central(2:end,1:2) - mask_central(1:end-1,1:2);
% delete outliers
diff_mask_central(abs(diff_mask_central(:,1))>34,1) = 0;
diff_mask_central(abs(diff_mask_central(:,2))>34,2) = 0;

% column 3 is the 2-norm shift distance considering both x,y direction
diff_mask_central(:,3) = sqrt(diff_mask_central(:,1).^2+diff_mask_central(:,2).^2);

%column 4 is the vdeio timestamp
diff_mask_central(:,4) = mask_central(2:end,3);

%% read stage moving information and save it in 'diff_mask_central(:,5,6)' according to time stamps in 'diff_mask_central(:,4)'
% rescale 'stage_time' when last stage motion is out of time scale
if stage_time(end)>  video_timestamp_time(end)
    % rescale stage_time to timestamp_time(end) and keep 4 digits after point
    stage_time_adj = round(stage_time/(stage_time(end)+0.1)*video_timestamp_time(end),3);
else
    stage_time_adj = stage_time;
end

% insert the stage motion according to time
pt = 2;
for ii = 1:size(diff_mask_central,1);
    disp(ii)
    % insert the stage motion in the right time
    if (pt <= numel(stage_time)) && (stage_time_adj(pt)<diff_mask_central(ii,4))
        % estimate the stage motion in cvs
        diff_mask_central(ii,5:6) = stage_xy(pt,:)- stage_xy(pt-1,:);
        % find next shift
        pt = pt + 1;
    end
end
num_pt = pt;

%% normalize shift distance

% normalize standard stage moving
% please be care about the y_pixel and x_pixle here!!!
diff_mask_central(:,5) = diff_mask_central(:,5)/x_pixel_per_microns;
diff_mask_central(:,6) = diff_mask_central(:,6)/y_pixel_per_microns;
diff_mask_central(:,7) = sqrt(diff_mask_central(:,5).^2+diff_mask_central(:,6).^2);

if size(video_timestamp_ind,1)-1~=size(diff_mask_central,1)
    error('time stamp frame number error');
end
diff_mask_central(:,8) = video_timestamp_ind(2:end);

% fix the missing frames as 0
diff_mask_central_full = zeros(video_timestamp_ind(end),size(diff_mask_central,2));
exist_ind = video_timestamp_ind(2:end);
curr_ind1 = 1;
for ii = 1: size(diff_mask_central_full,1);
    if ii == exist_ind(curr_ind1)
        diff_mask_central_full(ii,:) = diff_mask_central(curr_ind1,:);
        curr_ind1 = curr_ind1 +1;
    else
        diff_mask_central_full(ii,:) = zeros(1,size(diff_mask_central_full,2));
        diff_mask_central_full(ii,1:4) = 1e-5;
    end
end


%{
% tell if any motion in csv is too large
if frame_total ~= length(real_time_frame)
    csv_large_ind = find(diff_mask_central(:,7)>80);
    error('some indexes in csv has larger magnitude, eg: \n %d=%d, etc.',csv_large_ind(1), diff_mask_central(csv_large_ind(1),7));
end
%}

%% build a match/alignment between peaks of csv and WormShift
diff_leng = size(diff_mask_central,1);
moving_frame = zeros(diff_leng,4);

% this if is VERY IMPORTANT, you can adjust threshold parameters here
% if it is swimming video, moving frames are only the ones with
% large center shift

if is_swimming>0
    moving_frame(:,1) = (diff_mask_central(:,3)>4.5);
else
    % if it is an 'on food' video, moving frames are ones with
    % large x/yShift, and good center shift value, or Very large
    % center shift value (too large, so it was set to 0 in previous codes)
    
    % (xShift|yShift>4) & (diff_central is large | diff_central is set to 0)
    moving_frame(:,1) = (((abs(xShift)>4)|(abs(yShift)>4))>0|(abs_diff_fra_trim>1))&...
        ((diff_mask_central(:,3)>1.8)|(diff_mask_central(:,1)==0)|(diff_mask_central(:,2)==0));
end
%moving_frame(:,2) = moving_frame(:,1);
moving_frame(:,2) = 0;

% calculate the stage motion length
motion_len = zeros(10,1);
% parameter to record the current motion stage lengh
curr_motion_len = 0;
% go through all indexes to find the average motion length
for qq = 2: diff_leng;
    if moving_frame(qq, 1) == 1
        curr_motion_len = curr_motion_len +1;
    elseif moving_frame(qq, 1) == 0 && moving_frame(qq-1, 1) == 1
        bot = min(10,curr_motion_len);
        motion_len(bot) = motion_len(bot)+1;
        curr_motion_len = 0;
    end
end
ave_motion_len = round(motion_len'*(1:10)'/sum( motion_len));
% reduce each peak to "single frame length"
for qq = 1: diff_leng-2;
    %                 % debug use
    %                 if qq == 42841
    %                     qq
    %                 end
    bot = max(qq-ave_motion_len,1);
    top = max(qq-1,1);
    if moving_frame(qq, 1) == 1 && sum(moving_frame(bot:top,2)) ==0 && moving_frame(qq, 2) == 0
        if all(moving_frame(qq+1:qq+2,1)==[0;0])
            moving_frame(qq, 2) = 1;
            %                     elseif moving_frame(qq+1:qq+3,1)==[1;1;1]
            %                         moving_frame(qq, 2) = 0;
        else
            moving_frame(qq:qq+1, 2)= [0;1];
        end
    end
end

% 'mov_fra_ind' is the potential moving frames
mov_fra_ind = (find(moving_frame(:,2)>0));
%mov_fra_ind_full = diff_mask_central((find(moving_frame(:,2)>0)),8);
% 'csv_ind' is the real moving frames imported from csv
csv_ind = (find(diff_mask_central(:,7)>0));
%csv_ind_full = diff_mask_central((find(diff_mask_central(:,7)>0)),8);

%% adjust the scale of "mov_fra_ind" to make it match to the "csv_ind"

match_mtx = zeros(20,40);
min_match_mtx_elem = 1e10;
min_fra_ind_match = [];

% mm1 varies from 1 to 20 peaks shift
for mm1 = 1:min(20,length(mov_fra_ind));
    % mm2 varies from 0.94 to 1.06, with 0.03 step size
    for mm2 = 1:40;% 20:20 
        % the key function to do the match job
        [match_mtx(mm1, mm2),mov_fra_ind_match] = cal_match_score(csv_ind, mov_fra_ind, mm1, 0.94+0.003*mm2, diff_mask_central(:,7), diff_mask_central(:,3));

        
        if match_mtx(mm1, mm2) < min_match_mtx_elem
            min_match_mtx_elem = match_mtx(mm1, mm2);
            min_fra_ind_match = mov_fra_ind_match;
            % save the best shift distance and scale parameter
            % mm1 is the shift distance
            mm1_min_ind = mm1;
            % mm2 is the scale parameter
            mm2_min_ind = mm2;
        end
    end
end
% calculate the orignal indexes, and sort it to ascending order
min_fra_ind_match(:,5) = sort(round((min_fra_ind_match(:,2)-csv_ind(1))/(0.94+0.003*mm2_min_ind)+mov_fra_ind(mm1_min_ind)));
min_fra_ind_match(:,6) = sort(round((min_fra_ind_match(:,1)-csv_ind(1))/(0.94+0.003*mm2_min_ind)+mov_fra_ind(mm1_min_ind)));

%% refine aligment indexes

% gaussian window covers the impluse peak
% gaussian window size
filter_length = 71;
% half of the window size
half_filter = (filter_length-1)/2 ;
% generate the gausssian window vector


%code of gausswin, otherwise one requires the signal processing toolbox. guass_window = gausswin(filter_length);
L = filter_length; a = 2.5;
N = L-1;
n = (0:N)'-N/2;
guass_window = exp(-(1/2)*(a*n/(N/2)).^2);

% parameters to calculate the 'shift_to_left'
shift_weights = exp(-2.4:0.6:0)';
normal_shift_weights = shift_weights/sum(shift_weights);
% only conder nearest 10 peaks to calculate the 'shift_to_left'
shift_consider = zeros(1,length(shift_weights));

mov_compare = zeros(diff_leng,2);
len_mov_fra_ind_extend = length(mov_fra_ind);
% estimate the shift distance around each frame, +-3 frames
for ll = 1:len_mov_fra_ind_extend
    mov_compare(mov_fra_ind(ll),1:2) = sum(diff_mask_central(max(1,mov_fra_ind(ll)- ave_motion_len):min(diff_leng,mov_fra_ind(ll)+ ave_motion_len),1:2));
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
    if ii == 51
        ii
    end
    
    %ind_consider0 = round((min_fra_ind_match(ii,1)+min_fra_ind_match(ii,5))/2)-shift_to_left(ii);
    ind_consider0 = round(min_fra_ind_match(ii,6))-shift_to_left(ii);
    
    % actually it should be: start_in_abs_diff = ind_consider -
    % half_filter+ half_filter; to find the indexes of values of interests
    % in 'abs_diff'
    start_in_abs_diff0 = max(1,ind_consider0-half_filter);
    % consider the case when 'end_inabs_diff0' is larger than 'diff_leng'
    end_in_abs_diff0 = min(diff_leng,ind_consider0+half_filter);
    if start_in_abs_diff0 > end_in_abs_diff0
        sum_abs_diff0_sec = sum(end_in_abs_diff0 :end_in_abs_diff0+half_filter);
        ind_in_consider = end_in_abs_diff0:end_in_abs_diff0+half_filter;
    else
        % find the peak index with maximum weights, within the areas of 'positive/negtive areas are both large'
        sum_abs_diff0_sec = sum(ind_consider0:end_in_abs_diff0+half_filter);
        ind_in_consider = ind_consider0:end_in_abs_diff0+half_filter;
    end
    
    
    % use a guassian window to esitmate the best weights
    peak_weights0 = ones(length(ind_in_consider),1)*diff_mask_central(min_fra_ind_match(ii,1),5:6)-diff0(ind_in_consider,1:2).*(guass_window(round((filter_length-length(ind_in_consider))/2)+(1:length(ind_in_consider)))*ones(1,2));
    % debug use
    if sum_abs_diff0_sec==0
        start_in_abs_diff0
        end_in_abs_diff0
        error('no impluses in the interval')
    end
    [weigths_max0, weights_max_ind0] = min(sum(abs(peak_weights0),2));
    %  results of indexes and absolute values of diff
    mask_ind(ii,1) = ind_in_consider(weights_max_ind0)-half_filter ;
    mask_ind(ii,2) = sum(diff_mask_central(max(1,mask_ind(ii,1)-ave_motion_len):min(diff_leng,mask_ind(ii,1)+ave_motion_len),3));
    
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
% shift weights in choosing the penalty of corresponding peaks
shift_para = 100;
cross_para = 0.1;
min_fra_ind_match(:,7) = zeros(size(min_fra_ind_match(:,6)));
for iin = 1: no_nonzero_csv_ind;
    if iin == 95
        iin
    end
    % compare two alignment results, if they are the same, then choose any
    % one of them
    if (min_fra_ind_match(iin,5) ==  min_fra_ind_match(iin,6)) || (iin == 1)
        min_fra_ind_match(iin,7) = min_fra_ind_match(iin,6);
    else
        % if two alignments are not the same, compare the shift in x,y axis,
        % in terms of 'shift of central of area' and 'result of cross
        % correlation', respectively.
        
        % magnitudes of shift of area centre
        xy_sum_diffCentr_5 = sum(diff_mask_central(max(1,min_fra_ind_match(iin,5)- ave_motion_len):...
            min(diff_leng,min_fra_ind_match(iin,5)+ ave_motion_len),1:2),1)-...
            (diff_mask_central(max(1,min_fra_ind_match(iin,5)- ave_motion_len-1),1:2)+...
            diff_mask_central(min(diff_leng,min_fra_ind_match(iin,5)+ ave_motion_len+1),1:2));
        % magnitudes of shift of area centre
        xy_sum_diffCentr_6 = sum(diff_mask_central(max(1,min_fra_ind_match(iin,6)- ave_motion_len):...
            min(diff_leng,min_fra_ind_match(iin,6)+ ave_motion_len),1:2),1)-...
            (diff_mask_central(max(1,min_fra_ind_match(iin,6)- ave_motion_len-1),1:2)+...
            diff_mask_central(min(diff_leng,min_fra_ind_match(iin,6)+ ave_motion_len+1),1:2));
        % magnitudes of the sum of crossCorrelation around each peak, for
        % the 5th column of 'min_fra_ind_match'
        xy_sum_CrossCor_5 = [sum(xShift(max(1,min_fra_ind_match(iin,5)- ave_motion_len):...
            min(diff_leng,min_fra_ind_match(iin,5)+ ave_motion_len))),...
            sum(yShift(max(1,min_fra_ind_match(iin,5)- ave_motion_len):...
            min(diff_leng,min_fra_ind_match(iin,5)+ ave_motion_len)))];
        % magnitudes of the sum of crossCorrelation around each peak, for the 6th column of 'min_fra_ind_match'
        xy_sum_CrossCor_6 = [sum(xShift(max(1,min_fra_ind_match(iin,6)- ave_motion_len):...
            min(diff_leng,min_fra_ind_match(iin,6)+ ave_motion_len))),...
            sum(yShift(max(1,min_fra_ind_match(iin,6)- ave_motion_len):...
            min(diff_leng,min_fra_ind_match(iin,6)+ ave_motion_len)))];
        xy_csv = diff_mask_central(min_fra_ind_match(iin,1),5:6);
        % estimate the error in terms of col
        error_5 = sum(abs(xy_sum_diffCentr_5 - xy_csv)+cross_para*abs(xy_sum_CrossCor_5 - xy_csv))...
            + shift_para*abs((min_fra_ind_match(iin,5)-min_fra_ind_match(iin-1,7))/...
            (min_fra_ind_match(iin,1)-min_fra_ind_match(iin-1,1))-1);
        error_6 = sum(abs(xy_sum_diffCentr_6 - xy_csv)+cross_para*abs(xy_sum_CrossCor_6 - xy_csv))...
            + shift_para*abs((min_fra_ind_match(iin,6)-min_fra_ind_match(iin-1,7))/...
            (min_fra_ind_match(iin,1)-min_fra_ind_match(iin-1,1))-1);
        
        % select the alignment result with smaller x,y shift errors and
        % non-repeated indexes.
        if (error_5 <= error_6) && ~any(abs(min_fra_ind_match(1:iin,7)-min_fra_ind_match(iin,5))<1e-5)
            min_fra_ind_match(iin,7) = min_fra_ind_match(iin,5);
        elseif (error_6 < error_5) && ~any(abs(min_fra_ind_match(1:iin,7)-min_fra_ind_match(iin,6))<1e-5)
            min_fra_ind_match(iin,7) = min_fra_ind_match(iin,6);
        else
            back_ind = 1;
            min_fra_ind_match(iin,7) = min_fra_ind_match(iin,5);
            % if there is a duplicated index, choose the 5th column as the
            % result
            while length(unique(min_fra_ind_match(1:iin,7)))~=iin && iin-back_ind > 0
                min_fra_ind_match(iin-back_ind,7)=min_fra_ind_match(iin-back_ind,5);
                % the backwards index, increase from 1, until no duplicated
                % index exist in final index result vector
                back_ind = back_ind+1;
            end
        end
    end
    
end
min_fra_ind_match_7_diff0 = min_fra_ind_match(:,7)-min_fra_ind_match(:,1);
min_fra_ind_match_7_diff = min_fra_ind_match_7_diff0(2:end)-min_fra_ind_match_7_diff0(1:end-1);
if sum(abs(abs(min_fra_ind_match_7_diff)-mean(min_fra_ind_match_7_diff))>50)>=2
    mm1_min_ind
    mm2_min_ind
    min_fra_ind_match(:,7) = min_fra_ind_match(:,5);
end

%% important: obtain the alighment indexes here

stage_move_x = zeros(diff_leng,1);
stage_move_y = zeros(diff_leng,1);

range1 = 15;

% align csv_ind to mask_ind
for ii = 1:no_nonzero_csv_ind;
    % calculate stage moving vectors
    stage_move_x = align_func(stage_move_x, diff_mask_central(:,1), diff_mask_central(:,5),range1,  min_fra_ind_match(:,7), csv_ind, ii);
end
for ii = 1:no_nonzero_csv_ind;
    % calculate stage moving vectors
    stage_move_y = align_func(stage_move_y, diff_mask_central(:,2), diff_mask_central(:,6),range1,  min_fra_ind_match(:,7), csv_ind, ii);
end

diff_mask_central(:,9) = stage_move_x;
diff_mask_central(:,10) = stage_move_y;
% end


%% show all skeletons
% location of stage
loc_stage = zeros(size(diff_mask_central,1),2);
loc_stage(min_fra_ind_match(:,7),1:2) = diff_mask_central(min_fra_ind_match(:,1),5:6);

stage_mov_x_cum = cumsum(loc_stage(:,1));
stage_mov_y_cum = cumsum(loc_stage(:,2));

% tell the length of each peak/stage_motion
moving_frame(min_fra_ind_match(:,7),3)=1;
sese = ones(ave_motion_len+1,1);

% dilate each peak to a length of 'ave_motion_len+1'
moving_frame(:,3)=imdilate(moving_frame(:,3),sese);  % also can consider '(conv(moving_frame(:,3),sese))>0'

% adjust the threshold of abs_diff_fra to have a tighter peak length
abs_diff_fra_thres2 = abs_diff_fra;
abs_diff_fra_thres2(abs_diff_fra<(graythresh(abs_diff_fra/max(abs_diff_fra))*max(abs_diff_fra)*0.7))=0;  % no 0.9 compensate here

% generate a new 'moving_frame(:,1)'
%moving_frame1 =  (((abs(xShift)>4)|(abs(yShift)>4))>0|(abs_diff_fra_thres2>1));
%moving_frame1 =  (abs_diff_fra_thres2>1);
moving_frame1 = (((abs(xShift)>4)|(abs(yShift)>4))>0|(abs_diff_fra_thres2>1))&...
        ((diff_mask_central(:,3)>2)|(diff_mask_central(:,1)==0)|(diff_mask_central(:,2)==0));
    
moving_frame(:,4) = moving_frame1.*moving_frame(:,3);
last_peak = find(moving_frame(:,3)==1, 1, 'last');%max(find(moving_frame(:,3)==1));
moving_frame(min(last_peak,size(moving_frame,1)):end,4)= moving_frame(min(last_peak,size(moving_frame,1)):end,1);

% fix small problems in 'moving_frame(:,4)'
for repeat_i =1:3;   % repeat this process twice
    temp_moving_fra = moving_frame1 + moving_frame(:,4);
    for moving_ii = 2:diff_leng-4;
        if all(temp_moving_fra(moving_ii:moving_ii+1)==[1;2]) && (moving_frame1(moving_ii)==1) % extend length
            moving_frame(moving_ii:moving_ii+1,4) = [1;1];
        elseif all(temp_moving_fra(moving_ii:moving_ii+1) == [2;1]) && (moving_frame1(moving_ii+1)==1)
            moving_frame(moving_ii:moving_ii+1,4) = [1;1];
        elseif all(temp_moving_fra(moving_ii:moving_ii+2) == [2;0;2])  % fix a gap in length
            moving_frame(moving_ii:moving_ii+2,4) = [1;1;1];
        end
    end
end
%%
cancel_fra_ind2 = find(moving_frame(:,4) ==1);
cancel_fra_ind = cancel_fra_ind2+1;

skeleton_hdf5 = h5read(skeletons_file,'/skeleton');
x_ske = (reshape(skeleton_hdf5(1,:,:), size(skeleton_hdf5,2),size(skeleton_hdf5,3)))';
y_ske = (reshape(skeleton_hdf5(2,:,:), size(skeleton_hdf5,2),size(skeleton_hdf5,3)))';
x_ske_cum = x_ske;
y_ske_cum = y_ske;
    x_ske_cum(2:end,:) = x_ske(2:end,:) - stage_mov_x_cum*ones(1,size(x_ske,2));
    y_ske_cum(2:end,:) = y_ske(2:end,:) - stage_mov_y_cum*ones(1,size(y_ske,2));
    cancel_fra_ind = cancel_fra_ind2+1;
    x_ske_cum(cancel_fra_ind,:) = NaN;
    y_ske_cum(cancel_fra_ind,:) = NaN;
left_fra_ind = setdiff([1:size(skeleton_hdf5,2)], cancel_fra_ind);
    

stage_motion_x = [stage_mov_x_cum(1);stage_mov_x_cum];
stage_motion_y = [stage_mov_y_cum(1);stage_mov_y_cum];
stage_motion_x(cancel_fra_ind) = NaN;
stage_motion_y(cancel_fra_ind) = NaN;
stage_vec = [stage_motion_x';stage_motion_y'];

pixel_to_micro = [-(x_pixel_per_microns), -(y_pixel_per_microns)];

%%save stage vector
fid = H5F.open(skeletons_file,'H5F_ACC_RDWR','H5P_DEFAULT');
if H5L.exists(fid,'stage_vec','H5P_DEFAULT')
    H5L.delete(fid,'stage_vec','H5P_DEFAULT');
end
H5F.close(fid);
h5create(skeletons_file, '/stage_vec', size(stage_vec), 'Datatype', 'double', ...
'Chunksize', size(stage_vec), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
h5write(skeletons_file, '/stage_vec', stage_vec);


pixels2microns_x = h5readatt(masked_image_file, '/mask', 'pixels2microns_x');
pixels2microns_y = h5readatt(masked_image_file, '/mask', 'pixels2microns_y');
h5writeatt(skeletons_file, '/trajectories_data', 'pixels2microns_x', pixels2microns_x);
h5writeatt(skeletons_file, '/trajectories_data', 'pixels2microns_y', pixels2microns_y);

%% save the stage motion vector
% stage_vec_save = [trajectories_folder,name,'_stage_vec.mat'];
% save(stage_vec_save,'stage_vec');

%% compare to the 'features.mat' (optional)
 % CompareToFeature_func_ver2

%%
% function CompareToFeature_func
% features_mat = [real_features_folder,name_temp,'_features.mat'];
load(features_mat);
% N:\Andre\results-12-05-10\Laura Grundy\egl-17\e1313\CB1313\on_food\XX\30m_wait\L\tracker_2\2010-07-09___11_43_13\

% index of not-NaN entries in '.mat' x coordinates
Nannan_ind = worm.posture.skeleton.x(1,:)>0;
Nan2_ind = ~(info.video.annotations.frames==3);

worm_skeleton_x = worm.posture.skeleton.x(:,Nannan_ind);
worm_skeleton_y = worm.posture.skeleton.y(:,Nannan_ind);

worm_skeleton_x2 = worm.posture.skeleton.x(:,Nan2_ind);
worm_skeleton_y2 = worm.posture.skeleton.y(:,Nan2_ind);

tt_step = 15;

%%%%%%
%left_fra_sta_ind = diff_mask_central(cancel_fra_ind2,8)+1;
%left_fra_sta_ind = diff_mask_central(left_fra_ind-1,8)+1;

x_ske_cum_left = x_ske_cum(left_fra_ind,:)*(-x_pixel_per_microns);  % x_ske_cum
y_ske_cum_left = y_ske_cum(left_fra_ind,:)*(-y_pixel_per_microns);  % y_ske_cum

Nannan_ind2 = ~isnan(x_ske_cum_left(:,1));
% adjust worm skeleton based on the first non-nan entry in
% worm_skeleton_x,y
for ii = 1:1e3;
    if ~isnan(x_ske_cum_left(ii,1))
        xy_adjust_ind = ii;
        break;
    end
end
x_ske_cum_left2 = x_ske_cum_left(Nannan_ind2,:)+(worm_skeleton_x(xy_adjust_ind,1)-x_ske_cum_left(xy_adjust_ind,1));
y_ske_cum_left2 = y_ske_cum_left(Nannan_ind2,:)+(worm_skeleton_y(xy_adjust_ind,1)-y_ske_cum_left(xy_adjust_ind,1));


fig66 = figure(66),plot(worm_skeleton_x(:,1),worm_skeleton_y(:,1),'r' ); axis equal; hold on

hold on
%fig79 = figure(79),

for tt_1 = 1:tt_step:size(worm_skeleton_x,2);
    %if ismember(tt_1+1, cancel_fra_ind)
    if ~(isnan(sum(worm_skeleton_x(:,tt_1))) | isnan(sum(worm_skeleton_y(:,tt_1))))
        plot(worm_skeleton_x(:,tt_1),worm_skeleton_y(:,tt_1),'r' );
    end
end

plot(x_ske_cum_left2(1,:),y_ske_cum_left2(1,:),'b' ); axis equal; 

for tt_1 = 1:tt_step:size(x_ske_cum_left2,1);
    plot(x_ske_cum_left2(tt_1,:),y_ske_cum_left2(tt_1,:) ,'b');
end
hold off,

frame_total = size(video_timestamp_ind,1);
frame_total_ske = frame_total-sum(isnan(x_ske_cum_left(:,1)));
% show frames shown
frame_show1 = size(x_ske_cum_left2,1);

%%%%%%%%%%%%%%

%fig67 = figure(67),plot(worm_skeleton_x(:,size(worm_skeleton_x,2)-200),worm_skeleton_y(:,size(worm_skeleton_x,2)-200) ); axis equal; hold on

% show frames shown
fra_total_ori = size(worm.posture.skeleton.x,2)-sum(info.video.annotations.frames == 3);
fra_total_seg_ori = sum(info.video.annotations.frames == 1)+(sum(info.video.annotations.frames == 2));
frame_show_ori = size(worm_skeleton_x,2);
%%%%%

%savefig(fig66,[trajectories_folder,name,'-align_ske.fig']);
% savefig(fig67,[trajectories_folder,name,'-real_ske.fig']);

% show several frame number
frame_total     % existing frames (all - missing frames)
frame_total_ske % frames with skeleton (all - missing frames - cannot recognize skeleton)
frame_show1     % shown frames (all - missing frames - cannot recognize skeleton - stage motion )
fra_total_ori   % total original (all - reference 3)
fra_total_seg_ori  %   (ref 1 + ref 2)
frame_show_ori     %    (ref 1 )


% calculate the skeleton difference
x_ske_full = x_ske_cum*(-x_pixel_per_microns)+(worm_skeleton_x2(xy_adjust_ind,1)-x_ske_cum_left(xy_adjust_ind,1));
y_ske_full = y_ske_cum*(-y_pixel_per_microns)+(worm_skeleton_y2(xy_adjust_ind,1)-y_ske_cum_left(xy_adjust_ind,1));

min_full_ind = min(fra_total_ori,frame_total-1);

diff_x_row = sum(abs(x_ske_full(1:min_full_ind,:)-(worm_skeleton_x2(:,1:min_full_ind))'),2);
diff_y_row = sum(abs(y_ske_full(1:min_full_ind,:)-(worm_skeleton_y2(:,1:min_full_ind))'),2);
sum_diff = diff_x_row+diff_y_row;
%sum_diff = diff_x_row(~isnan(diff_x_row))+diff_y_row(~isnan(diff_y_row));
figure(55), plot(sum_diff)
title('Absolute difference of skeleton between results and segworm ');

centr_diff = sum(abs(x_ske_full(2:end,:)-x_ske_full(1:end-1,:))+abs(y_ske_full(2:end,:)-y_ske_full(1:end-1,:)),2);
centr_diff_seg = sum(abs(worm_skeleton_x2(:,2:end)-worm_skeleton_x2(:,1:end-1))+abs(worm_skeleton_y2(:,2:end)-worm_skeleton_y2(:,1:end-1)),1);
figure(56),  plot(centr_diff,'b'), hold on,plot(centr_diff_seg,'r'), hold off
title('Comparison of skeleton differences of results and segworm along time, respectively ');

% max_gap_shift = max(abs(gap_shift));
% txt_name = [real_features_folder,name_temp, '-summary.txt'];
% fid_txt = fopen(txt_name,'wt');
% fprintf(fid_txt, 'max_gap_shift = %d \n gap_shift frame number of results: \n frame_total =%d \n frame_total_ske = %d \n frame_show =%d \n fra_total_ori =%d \n fra_total_seg_ori = %d \n frame_show_ori =%d \n skeleton difference: %d',...
%     max_gap_shift, frame_total_ske,frame_total_ske, frame_show1,fra_total_ori, fra_total_seg_ori,frame_show_ori, sum_diff);
% fclose(fid_txt);


figure(71), plot(info.video.annotations.frames == 2,'r')
peak_ind_in_total = zeros(length(info.video.annotations.frames),1); 
peak_ind_in_total(diff_mask_central(cancel_fra_ind2,8)+1) = 1.15;
%hold on, plot(peak_ind_in_total,'r')
hold on, plot(peak_ind_in_total,'b')
hold off,
