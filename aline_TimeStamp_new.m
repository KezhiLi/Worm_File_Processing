% Align the stage moving info in .csv file, and video file
% 
% 
% 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 10/12/2015
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

% 
% %% file root preparation
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
% 
% diff_mask_central = mask_central(2:end,1:2) - mask_central(1:end-1,1:2);
% diff_mask_central(:,3) = sqrt(diff_mask_central(:,1).^2+diff_mask_central(:,2).^2);
% diff_mask_central(:,4) = mask_central(2:end,3);
% 
% %% read stage moving information and save it in 'diff_mask_central(:,5,6)' according to time stamps in 'diff_mask_central(:,4)'
% pt = 2;
% for ii = 1:frame_total-1;
%     ii
%     if (pt <= data_rows) & (media_time_vec(pt)<diff_mask_central(ii,4))
%           % switch x,y axis 
%         diff_mask_central(ii,6:-1:5) = stage_xy(pt,:)- stage_xy(pt-1,:);
%         pt = pt + 1;
%     end 
% end
% % calculate the absolute moving distance based on x-axis and y-axis
% diff_mask_central(:,7) = sqrt(diff_mask_central(:,5).^2+diff_mask_central(:,6).^2);
% 
% 
% % %% show some align figures
% % figure, plot(diff_mask_central(1:10000,3)), hold on, plot(diff_mask_central(1:10000,7)/25,'red')
% % figure, plot(diff_mask_central(10000:20000,3)), hold on, plot(diff_mask_central(10000:20000,7)/25,'red')
% % figure, plot(diff_mask_central(end-10000:end,3)), hold on, plot(diff_mask_central(end-10000:end,7)/25,'red')
% 
%  save example1.mat
% 
load example1.mat

%% read pixels per microns from xml files
for qq = 1:2;
    if qq == 1
        findLabel1 = 'microns';
    elseif qq == 2
        findLabel1 = 'pixels';
    end

    % end - 8 because: excel name ends with '.log.csv', 8 letters
    xDoc1 = xmlread([folder,excel_name(1:end-8),'.info.xml']);
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
diff_mask_central(:,8) = pixel_diff(:,3)/10e9;
sorted_diff_area = sort(diff_mask_central(:,8),'descend');
% maximum keep 15% time slots that are potentially belongs to stage moving
thresh_diff_area = sorted_diff_area(round(0.15*size(sorted_diff_area,1)));

%% find corresponding impulses, for x
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
abs_diff = abs([zeros(half_filter,1);diff_mask_central(:,1);zeros(half_filter,1)]);
% possible absolute difference based on "difference multiplication"
possible_abs_diff =  [zeros(half_filter,1);(diff_mask_central(:,8)>thresh_diff_area);zeros(half_filter,1)];
% the indexes of nonzero peaks in csv
csv_ind_x =find( abs(diff_mask_central(:,5))>0);
% number of nonzeros peaks
no_nonzero = size(csv_ind_x,1);

diff_mask_central(:,9:10) = zeros(size(diff_mask_central,1),2);

% determine x mask index
for ii = 1:no_nonzero;
    % abstract 'shift_to_left' to find the most possible csv index
    ind_consider = csv_ind_x(ii)-shift_to_left;
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
    shift_consider(1:end-1) = shift_consider(2:end);
    shift_consider(end) = csv_ind_x(ii)-mask_ind_x(ii,1);
    shift_to_left = round(shift_consider*normal_shift_weights);
    
end

%% find corresponding impulses, for y

% only conder nearest 10 peaks to calculate the 'shift_to_left'
shift_consider = zeros(1,10);

% compensate the difference: shift impulse index to left 
shift_to_left = 0;

% absolute difference of stage central with extensions on both sides with
% all 0. 
abs_diff = abs([zeros(half_filter,1);diff_mask_central(:,2);zeros(half_filter,1)]);
% the indexes of nonzero peaks in csv
csv_ind_y = find( abs(diff_mask_central(:,6))>0);
% number of nonzeros peaks
no_nonzero = size(csv_ind_y,1);

% determine y mask index
for ii = 1:no_nonzero;
    % abstract 'shift_to_left' to find the most possible csv index
    ind_consider = csv_ind_y(ii)-shift_to_left;
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
    shift_consider(end) = csv_ind_y(ii)-mask_ind_y(ii,1);
    shift_to_left = round(shift_consider*normal_shift_weights);
    
end


figure, plot(diff_mask_central(1:10000,3)), hold on, plot(diff_mask_central(1:10000,7)/25,'red')
%hold on, plot(mask_ind(1:10000,1),mask_ind(1:10000,2),'m*');
hold on, plot(pixel_diff(1:10000,3)/10e9,'go-')


%% calculate stage moving vectors
% range1 is the left and right range of searching area, in alignment. So
% the the window size will be 2*range1+1
range1 = 10;

% initilaize results vector
stage_move_x = zeros(size(diff_mask_central,1),1);
stage_move_y = zeros(size(diff_mask_central,1),1);
% align csv_ind to mask_ind
for ii = 1:no_nonzero;
    % calculate stage moving vectors
    stage_move_x = align_func(stage_move_x, diff_mask_central(:,1), diff_mask_central(:,5),range1, mask_ind_x, csv_ind_x, ii);
    stage_move_y = align_func(stage_move_y, diff_mask_central(:,2), diff_mask_central(:,6),range1, mask_ind_y, csv_ind_y, ii);
end

diff_mask_central(:,9) = stage_move_x;
diff_mask_central(:,10) = stage_move_y;

