% for pp =1:20000;
%     pp
%     if pp > 3700
%         pp
%     end
%     figure(12),imshow(mask(:,:,pp));
% end


% set parameters
timeDiff = 1; % how many frames between aligned images?
dS = 1; % pixel downsampling factor (2 means half size)

% estimate transformation from one image frame to another
No_mask = size(mask, 3);
xShift = NaN(No_mask-timeDiff, 1);
yShift = NaN(No_mask-timeDiff, 1);
for ii = 1+timeDiff:No_mask
    disp(ii/No_mask)
    
    frame_bef = mask(1:dS:end, 1:dS:end, ii);
    frame_aft = mask(1:dS:end, 1:dS:end, ii - timeDiff);
    
    frame_sum = abs(frame_bef)+abs(frame_aft)>0;
    frm_sum_col = sum(frame_sum,1);
    frm_sum_row = sum(frame_sum,2);
    
    worm_ind_col = find(frm_sum_col>0);
    worm_ind_row = find(frm_sum_row>0);
    
    pix_buffer = 5;
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
    xShift(ii - timeDiff) = transMat.T(3, 1);
    yShift(ii - timeDiff) = transMat.T(3, 2);
    
end