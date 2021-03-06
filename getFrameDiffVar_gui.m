function  img_abs_diff = getFrameDiffVar_gui(masked_image_file, hObject, handles)

global terminated;
terminated = 0; 

%% mask information
mask_info = h5info(masked_image_file, '/mask');
% size of each frame
frame_size = mask_info.Dataspace.Size(1:2);
frame_total = mask_info.Dataspace.Size(3);

%% Calculate the variance between consecutive images

img_abs_diff = zeros(frame_total-1,1);

blk = 500;
frame_prev =  h5read(masked_image_file, '/mask', [1,1,1], [frame_size(1),frame_size(2), 1]);
cur_frame_blk = h5read(masked_image_file, '/mask', [1,1,1], ...
            [frame_size(1),frame_size(2), min(frame_total,blk)]);
for ii = 2:frame_total;
    if terminated ==1
        break;
    end
        
    % ii = nn*blk + mm
    nn = floor(ii/blk);
    mm = mod(ii, blk);
    mm_full = [mm, blk];
    
    % read hdf5 in blocks to save time, because h5read is time costy
    if mm == 1
        %progress
        fprintf('%2.2f%%\n', (ii+1)/frame_total*100)
        cur_frame_blk = h5read(masked_image_file, '/mask', [1,1,(nn*blk+mm)], ...
            [frame_size(1),frame_size(2), min(frame_total-nn*blk,blk)]);
        set(handles.text18,'string',['Frame Diff Variance progress: ', num2str((ii+1)/frame_total*100),'%']);
        drawnow()
    end
    
    
    frame_current = squeeze(cur_frame_blk(:,:,mm_full(1+(mm==0))));
    %frame_current = h5read(masked_image_file, '/mask', [1,1,ii], [frame_size(1),frame_size(2), 1]);
    
    % calculate the absolute difference between frames
    %img_abs_diff(ii-timeDiff) = sum(sum(abs(double(frame_bef) - double(frame_aft))));
    
    
    img_abs_diff(ii-1) = imMaskDiffVar(frame_prev,frame_current);
    %good = (frame_prev>0) & (frame_current>0);
    %A = var(double(frame_prev(good))-double(frame_current(good)), 1);
    %assert(abs(A - img_abs_diff(ii-1)) < 1e-2);
    
    %buffer to the previous frame
    frame_prev = frame_current;
end