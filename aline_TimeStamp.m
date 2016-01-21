% % csv file
% folder = 'N:\Kezhi\DataSet\AllFiles\VideosExl\nas207-1\from_pc207-7\copied_from_pc207-8\';
% excel_name ='240 CB4852 on food L_2011_02_24__17_14___3___11.log.csv'; 
% % hdf5 file 
% hdf_folder = 'N:\Kezhi\DataSet\AllFiles\MaskedVideos\nas207-1\from_pc207-7\copied_from_pc207-8\';
% hdf5_name = '240 CB4852 on food L_2011_02_24__17_14___3___11.hdf5';
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
% 
% diff_mask_central = mask_central(2:end,1:2) - mask_central(1:end-1,1:2);
% diff_mask_central(:,3) = sqrt(diff_mask_central(:,1).^2+diff_mask_central(:,2).^2);
% diff_mask_central(:,4) = mask_central(2:end,3);
% 
% pt = 2;
% for ii = 1:frame_total-1;
%     ii
%     if (pt <= data_rows) & (media_time_vec(pt)<diff_mask_central(ii,4))
%           % switch x,y axis 
%         diff_mask_central(ii,6:-1:5) = stage_xy(pt,:)- stage_xy(pt-1,:);
%         pt = pt + 1;
%     end 
% end
% % diff_mask_central(:,7) = sqrt(diff_mask_central(:,5).^2+diff_mask_central(:,6).^2);
% 
% 
% % show some align figures
% figure, plot(diff_mask_central(1:10000,3)), hold on, plot(diff_mask_central(1:10000,7)/25,'red')
% figure, plot(diff_mask_central(10000:20000,3)), hold on, plot(diff_mask_central(10000:20000,7)/25,'red')
% figure, plot(diff_mask_central(end-10000:end,3)), hold on, plot(diff_mask_central(end-10000:end,7)/25,'red')
% 
% save example1.mat
load example1.mat

%% normalize
x_pixel_per_microns = -99.29374099796074/20.9974;
y_pixel_per_microns = -100.31847305825985/20.9974;

% normalize standard stage moving
diff_mask_central(:,5) = diff_mask_central(:,5)/y_pixel_per_microns;
diff_mask_central(:,6) = diff_mask_central(:,6)/x_pixel_per_microns;

diff_mask_central(:,8) = pixel_diff(:,3)/10e9;

%% find corresponding impulses
% initialize important values
% compensate the difference: shift impulse index to left 
shift_to_left = 0;
% gaussian window covers the impluse peak
% gaussian window size
filter_length = 31;
% half of the window size
half_filter = (filter_length-1)/2 ;
% generate the gausssian window vector
guass_window = gausswin(filter_length);

% parameters to calculate the 'shift_to_left'
shift_weights = exp(-2.7:0.3:0)';
normal_shift_weights = shift_weights/sum(shift_weights);
% only conder nearest 10 peaks to calculate the 'shift_to_left'
shift_consider = zeros(1,10);

% absolute difference of stage central with extensions on both sides with
% all 0. 
abs_diff = [zeros(half_filter,1);diff_mask_central(:,3);zeros(half_filter,1)];
possible_abs_diff =  [zeros(half_filter,1);(diff_mask_central(:,8)>1);zeros(half_filter,1)];
% the indexes of nonzero peaks in csv
csv_ind =find( diff_mask_central(:,7)>0);
% number of nonzeros peaks
no_nonzero = size(csv_ind,1);

% results saved in mask_ind, 1st column is the index in csv, 2nd column is
% its absolute value of diff
mask_ind = size(no_nonzero,2);

for ii = 1:no_nonzero;
    % abstract 'shift_to_left' to find the most possible csv index
    ind_consider = csv_ind(ii)-shift_to_left;
    % actually it should be: start_in_abs_diff = ind_consider -
    % half_filter+ half_filter; to find the indexes of values of interests
    % in 'abs_diff'
    start_in_abs_diff = ind_consider;
    end_in_abs_diff = start_in_abs_diff+filter_length-1;
    % find the peak index with maximum weights, within the areas of 'positive/negtive areas are both large'
    peak_weights = abs_diff(start_in_abs_diff:end_in_abs_diff).*guass_window.*possible_abs_diff(start_in_abs_diff:end_in_abs_diff);
    [weigths_max, weights_max_ind] = max(peak_weights);
    % save results of indexes and absolute values of diff
    mask_ind(ii,1) = start_in_abs_diff  + weights_max_ind - 1 -half_filter;
    mask_ind(ii,2) = diff_mask_central(mask_ind(ii,1),3);
    % update 'shift_to_left'. it is caluclated based on last 10 shift
    % values
    shift_consider(1:end-1) = shift_consider(2:end);
    shift_consider(end) = csv_ind(ii)-mask_ind(ii,1);
    shift_to_left = round(shift_consider*normal_shift_weights);
    
end

figure, plot(diff_mask_central(1:10000,3)), hold on, plot(diff_mask_central(1:10000,7)/25,'red')
%hold on, plot(mask_ind(1:10000,1),mask_ind(1:10000,2),'m*');
hold on, plot(pixel_diff(1:10000,3)/10e9,'go-')

%% 
range1 = 10;
x_consider = zeros(2*range1+1,1);
x_align = zeros(2*range1+1,1);
% align csv_ind to mask_ind
for ii = 1:no_nonzero;
    pulse_cental_ii = mask_ind(ii,1);
    % consider -10 to 10, in total 21 points
    x_consider = diff_mask_central(pulse_cental_ii-range1:pulse_cental_ii+range1,1);
    % true stage moving pixels
    x_stagemove = diff_mask_central(csv_ind(ii),5);
    % set a threshold to distinguish impulse part and other part
    x_level = graythresh(x_consider);
    % find the central impulse indexes
    if x_consider(range1+1,1)>0
        impulse_ind_binary = x_consider>x_level;       
    else
        impulse_ind_binary = x_consider<x_level;    
    end
    % only keep the central impulses, cancel other fake ones, by
    % identifying the 2 closest change points near central point
    impulse_ind_binary_diff = abs(impulse_ind_binary(2:end) - impulse_ind_binary(1:end-1));
    if sum(impulse_ind_binary_diff)> 2.5
        [impulse_ind_binary_diff_val, impulse_ind_binary_diff_ind]=sort(impulse_ind_binary_diff,'descend');
        impluse_chg_ind = impulse_ind_binary_diff_ind(logical(impulse_ind_binary_diff_val));
        chg_pt1 = impluse_chg_ind(1);
        if chg_pt1 < range1+0.5
            chg_pt2 = impluse_chg_ind(impluse_chg_ind>range1+0.5);
            chg_pt2 =  chg_pt2(2);
        else
            chg_pt2 = chg_pt1;
            chg_pt1 = impluse_chg_ind(impluse_chg_ind<range1+0.5);
            chg_pt1 =  chg_pt1(2);
        end        
        impulse_ind_binary_new = zeros(2*range1,1);
        impulse_ind_binary_new(chg_pt1+1:chg_pt2) = 1;
    end
    x_align(impulse_ind_binary) = x_consider(impulse_ind_binary); 
    
end
