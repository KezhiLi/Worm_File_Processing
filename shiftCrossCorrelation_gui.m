function [xShift, yShift] = shiftCrossCorrelation_gui(masked_image_file, hObject, handles)
%
%
%
%
%
%
%

%% set parameters
timeDiff = 1; % how many frames between aligned images?
dS = 4; % pixel downsampling factor (2 means half size)

%% mask information
mask_info = h5info(masked_image_file, '/mask');
% size of each frame
frame_size = mask_info.Dataspace.Size(1:2);
frame_total = mask_info.Dataspace.Size(3);

%% estimate transformation from one image frame to another
xShift = NaN(frame_total-timeDiff, 1);
yShift = NaN(frame_total-timeDiff, 1);
% img_abs_diff = zeros(frame_total-timeDiff,1);

blk = 500;
frame_prev =  h5read(masked_image_file, '/mask', [1,1,1], [frame_size(1),frame_size(2), 1]);
cur_frame_blk = h5read(masked_image_file, '/mask', [1,1,1], ...
            [frame_size(1),frame_size(2), min(frame_total,blk)]);
for ii = 1+timeDiff:timeDiff:frame_total
    % ii = nn*blk + mm
    nn = floor(ii/blk);
    mm = mod(ii, blk);
    mm_full = [mm, blk];
    
    if mm == 1
        %progress
        fprintf('%2.2f%%\n', (ii+1)/frame_total*100)
        cur_frame_blk = h5read(masked_image_file, '/mask', [1,1,(nn*blk+mm)], ...
            [frame_size(1),frame_size(2),min(frame_total-nn*blk,blk)]);
        set(handles.text18,'string',['progress: ', num2str((ii+1)/frame_total*100),'%']);
        drawnow()
    end
    
    frame_current = squeeze(cur_frame_blk(:,:,mm_full(1+(mm==0))));
    %frame_current = h5read(masked_image_file, '/mask', [1,1,ii], [frame_size(1),frame_size(2), 1]);
    
    %downsampling 
    frame_bef  = frame_current(1:dS:end,1:dS:end);
    frame_aft = frame_prev(1:dS:end,1:dS:end);
    
    %buffer to the previous frame
    frame_prev = frame_current;
    
    
    % use sum to find the worms in 2 frame, and the background of value 0
    frame_sum = frame_bef + frame_aft > 0;
    frm_sum_col = sum(frame_sum,1);
    frm_sum_row = sum(frame_sum,2);
    
    % find the joint worm body
    worm_ind_col = find(frm_sum_col>0);
    worm_ind_row = find(frm_sum_row>0);
    
    if isempty(worm_ind_col) || isempty(worm_ind_row)
        %if there is an empty mask continue
        xShift(ii - timeDiff) = nan;
        yShift(ii - timeDiff) = nan;
        continue
    end
    
    % the pixel buffer round the worm body area
    pix_buffer = 5;
    
    % create a square of pixels that cover the worm body, only
    % focus on the change pixles to speed up the calculation
    
    bot_x = max(1,worm_ind_row(1)-pix_buffer);
    top_x = min(length(frm_sum_row) ,worm_ind_row(end)+pix_buffer);
    bot_y = max(1,worm_ind_col(1)-pix_buffer);
    top_y = min(length(frm_sum_col),worm_ind_col(end)+pix_buffer);
    
    frame_bef = frame_bef(bot_x:top_x, bot_y:top_y);
    frame_aft = frame_aft(bot_x:top_x, bot_y:top_y);
    
    % estimate shift between images
    transMat = imregcorr(frame_bef , frame_aft  , 'translation');
    
    % calculate the absolute difference between frames
    % img_abs_diff(ii-timeDiff) = sum(sum(abs( frame_bef - frame_aft)));
    
    % calculate the shift in x,y directions
    xShift(ii - timeDiff) = transMat.T(3, 1)*dS;
    yShift(ii - timeDiff) = transMat.T(3, 2)*dS;
    
end

% switch and reverse the xShift and yShift, due to the transform
% between space coodinates and matrix presentation
xShift_temp = xShift;
xShift = -yShift;
yShift = -xShift_temp;