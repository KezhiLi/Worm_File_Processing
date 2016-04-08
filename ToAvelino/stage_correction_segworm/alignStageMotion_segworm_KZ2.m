%function alignStageMotion(masked_image_file,skeletons_file, is_swimming)

%main_dir = '/Users/ajaver/Desktop/Videos/single_worm/agar_1/MaskedVideos/';
main_dir = 'Z:\single_worm\agar_2\MaskedVideos';
results_dir = strrep(main_dir, 'MaskedVideos', 'Results');
feat_dir = strrep(main_dir, 'MaskedVideos', 'Features');

% tell if it is a swimming video
is_swimming = false;

% problem:
% 1. one fake peak
% 2. one neglect peak
% 3. too many missing frames in somewhere in frame_diffs
% 4. too many noises(fake peaks) in frame_diffs (swimming video?)
% 5. single bad pixel globally
% 6. other problems: eg. fps error
% 7. segworm misses peaks
% 8. last extra peaks

files = dir(main_dir);
for iif = 1:numel(files);
    % bad files:
    % agar_2: 6(1),18(2),29(2,3,fail);
    % agar_1: 6(4,fail), 9(1), 26(5), 29(,)
    % agar_goa: 9(4,?), 10(6,fail), 24(7,fail), 26(8,..),
    
    file = files(iif);
    
    if ~isempty(regexp(file.name, '\w*.hdf5', 'ONCE'))
        
        clear is_stage_move movesI stage_locations
        
        fprintf('%i) %s\n', iif, file.name)
        masked_image_file = fullfile(main_dir, file.name);
        skeletons_file = fullfile(results_dir, strrep(file.name, '.hdf5', '_skeletons.hdf5'));
        features_mat = fullfile(feat_dir, strrep(file.name, '.hdf5', '_features.mat'));
        
        %% read time stamps. I should put this data into the masked files dir
        video_timestamp_ind = h5read(skeletons_file, '/timestamp/raw');
        video_timestamp_ind = video_timestamp_ind + 1; %correct for python indexing
        trajectories_data = h5read(skeletons_file, '/trajectories_data'); % trajectories data
        real_time_frame = trajectories_data.timestamp_time;
        mask_central = [trajectories_data.cnt_coord_x, trajectories_data.cnt_coord_y,real_time_frame];
        
        if any(isnan(video_timestamp_ind))
            disp('The timestamp is corrupt or do not exist')
            continue
        end
        
        video_timestamp_time = h5read(skeletons_file, '/timestamp/time');
        fps = 1/median(diff(video_timestamp_time));
        
        %% main matrix used in alignment
        diff_mask_central = zeros(size(mask_central,1)-1, 4);
        
        % column 1  and 2 is the difference of the central of the mask in x and y
        diff_mask_central(:, 1:2) = mask_central(2:end,1:2) - mask_central(1:end-1,1:2);
        % delete outliers
        diff_mask_central(isnan(diff_mask_central(:,1)),1) = 10;
        diff_mask_central(isnan(diff_mask_central(:,2)),2) = 10;
        
        % column 3 is the 2-norm shift distance considering both x,y direction
        diff_mask_central(:,3) = sqrt(diff_mask_central(:,1).^2+diff_mask_central(:,2).^2);
        diff_mask_central(:,4) = real_time_frame(2:end);
        
        %% Open the information file and read the tracking delay time.
        % (help from segworm findStageMovement)
        % 2. The info file contains the tracking delay. This delay represents the
        % minimum time between stage movements and, conversely, the maximum time it
        % takes for a stage movement to complete. If the delay is too small, the
        % stage movements become chaotic. We load the value for the delay.
        
        xml_info = h5read(masked_image_file, '/xml_info');
        %this is not the cleaneast but matlab does not have a xml parser from
        %text string
        dd = strsplit(xml_info, '<delay>');
        dd = strsplit(dd{2}, '</delay>');
        delay_str = dd{1};
        delay_time = str2double(delay_str) / 1000;
        delay_frames = ceil(delay_time * fps);
        
        %% Read the media times and locations from the log file.
        % (help from segworm findStageMovement)
        % 3. The log file contains the initial stage location at media time 0 as
        % well as the subsequent media times and locations per stage movement. Our
        % algorithm attempts to match the frame differences in the video (see step
        % 1) to the media times in this log file. Therefore, we load these media
        % times and stage locations.
        %from the .log.csv file
        stage_data = h5read(masked_image_file, '/stage_data');
        mediaTimes = stage_data.stage_time';%*60;
        locations = [stage_data.stage_x , stage_data.stage_y];
        
        %% Read the scale conversions, we would need this when we want to convert the pixels into microns
        pixelPerMicronX = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_x');
        pixelPerMicronY = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_y');
        
        normScale = sqrt((pixelPerMicronX ^ 2 + pixelPerMicronX ^ 2) / 2);
        pixelPerMicronScale =  normScale * [sign(pixelPerMicronX) sign(pixelPerMicronY)];
        
        % Compute the rotation matrix.
        %rotation = 1;
        angle = atan(pixelPerMicronY / pixelPerMicronX);
        if angle > 0
            angle = pi / 4 - angle;
        else
            angle = pi / 4 + angle;
        end
        cosAngle = cos(angle);
        sinAngle = sin(angle);
        rotation_matrix = [cosAngle, -sinAngle; sinAngle, cosAngle];
        
        %% calculate the variance of the difference between frames
        % Ev's code uses the full vectors without dropping frames
        % 1. video2Diff differentiates a video frame by frame and outputs the
        % differential variance. We load these frame differences.
        frame_diffs_d0 = getFrameDiffVar_fast(masked_image_file);
        
        
        %% KZ added
        
        % %calculate shift from cross correlation between frames, and get the absolute difference between images
        [xShift, yShift] = shiftCrossCorrelation_fast(masked_image_file);
        
        % check if the number of frames
        if size(diff_mask_central,1)~=length(xShift)
            disp('hdf5 and skeleton has different numbers of frames: hdf5: %d, skeleton: %d',...
                length(xShift)+1,size(diff_mask_central,1)+1);
            continue
        end
        % subsample rate in calculating xShift,yShift, used as a threshold
        % in 'hybrid_mask1'
        dS = 4;
        %       % logical opration mask to tune the frame_diff
        hybrid_mask0 = (diff_mask_central(:,3)>3)...%|(diff_mask_central(:,1)==0)|(diff_mask_central(:,2)==0)...
            |((abs(xShift)>dS)|(abs(yShift)>dS));
        % extend mask to include the edges of peaks
        hybrid_mask1 = imdilate(hybrid_mask0,[1;1;1]);
        % use frame_diffs_d
        frame_diffs_d = frame_diffs_d0;
        % estimate global threshold
        g_thre = 0.4*(graythresh(frame_diffs_d0/max(frame_diffs_d0))*max(frame_diffs_d0)+median(frame_diffs_d0));
        % estimate cancel indexes
        canc_ind = (frame_diffs_d >g_thre)&(hybrid_mask1<0.5);
        % cancel the peaks that are considered as fake peaks
        %frame_diffs_d(hybrid_mask1<0.5)=1e-5;
        frame_diffs_d(canc_ind) = frame_diffs_d(canc_ind)/5;
        % fix mistake in hybrid_mask in there is gap in side an interval
        for ii = 3:length(frame_diffs_d)-2;
            if (frame_diffs_d(ii)<1)&&(frame_diffs_d(ii-1)>1)&&(frame_diffs_d(ii+1)>1)
                frame_diffs_d(ii) = frame_diffs_d0(ii);
            elseif  frame_diffs_d(ii)>g_thre&(frame_diffs_d0(ii-2)<g_thre)&(frame_diffs_d0(ii-1)<g_thre)...
                    &(frame_diffs_d0(ii+1)<g_thre)&(frame_diffs_d0(ii+2)<g_thre)
                frame_diffs_d(ii) = frame_diffs_d(ii)/5;
            end
        end
        
        %% The shift makes everything a bit more complicated. I have to remove the first frame, before resizing the array considering the dropping frames.
        
        if numel(video_timestamp_ind) > numel(frame_diffs_d) + 1
            %i can tolerate one frame (two with respect to the frame_diff)
            %extra at the end of the timestamp
            video_timestamp_ind = video_timestamp_ind(1:numel(frame_diffs_d)+1);
        end
        
        frame_diffs = nan(1, max(video_timestamp_ind)-1);
        dd = video_timestamp_ind-min(video_timestamp_ind);
        dd = dd(dd>0);
        if numel(frame_diffs_d) ~= numel(dd)
            continue
        end
        frame_diffs(dd) = frame_diffs_d;
        
        %%
        stage_locations = [];
        
        % IMPORTANT: window_weights = [0~0.8]; the larger, the more weights given to
        % the central part of interval in consideration
        % eg. choose 0.1 for normal, 0.7 for extream case ;
        for  wind_weights = 0.1:0.2:0.9;
            try
                clear is_stage_move movesI stage_locations ME
                [is_stage_move, movesI, stage_locations] = findStageMovement_gs(frame_diffs, mediaTimes, locations, delay_frames,...
                    fps, wind_weights);
            catch ME
                fprintf('%i) %s\n', iif, file.name)
                disp(ME)
                is_stage_move = ones(1, numel(frame_diffs)+1);
                stage_locations = [];
                continue
            end
            if ~(exist('ME','var'))
                disp(['finally successful when wind_weights = ',num2str(wind_weights,'%.3f')])
                break;
            end
        end
        %%
        stage_vec = nan(numel(is_stage_move),2);
        if numel(movesI) == 2 && all(movesI==0)
            %there was no movements
            stage_vec(:,1) = stage_locations(1);
            stage_vec(:,2) = stage_locations(2);
            
        else
            %convert output into a vector that can be added to the skeletons file to obtain the real worm displacements
            
            for kk = 1:size(stage_locations,1)
                bot = max(1, movesI(kk,2)+1);
                top = min(numel(is_stage_move), movesI(kk+1,1)-1);
                stage_vec(bot:top, 1) = stage_locations(kk,1);
                stage_vec(bot:top, 2) = stage_locations(kk,2);
            end
        end
        
        %the nan values must match the spected video motions
        assert(all(isnan(stage_vec(:,1)) == is_stage_move))
        
        %prepare vectors to save into the hdf5 file.
        %Go back to the original movie indexing. I do not want to include the missing frames at this point.
        frame_diffs_d = frame_diffs_d';
        is_stage_move_d = int8(is_stage_move(video_timestamp_ind))';
        
        
        %% change into a format that i can add directly to the skeletons in skeletons_file
        stage_vec_d = stage_vec(video_timestamp_ind, :);
        
        %stage_vec_d(:,1) = stage_vec_d(:,1)*pixels2microns_y;
        %stage_vec_d(:,2) = stage_vec_d(:,2)*pixels2microns_x;
        stage_vec_d = stage_vec_d';
        
        %%
        %this removes crap from previous analysis
        %%save stage vector
        fid = H5F.open(skeletons_file,'H5F_ACC_RDWR','H5P_DEFAULT');
        if H5L.exists(fid,'/stage_vec2','H5P_DEFAULT')
            H5L.delete(fid,'/stage_vec2','H5P_DEFAULT');
        end
        
        if H5L.exists(fid,'/is_stage_move2','H5P_DEFAULT')
            H5L.delete(fid,'/is_stage_move2','H5P_DEFAULT');
        end
        H5F.close(fid);
        
        
        %% delete data from previous analysis if any
        fid = H5F.open(skeletons_file,'H5F_ACC_RDWR','H5P_DEFAULT');
        if H5L.exists(fid,'/stage_movement2','H5P_DEFAULT')
            gid = H5G.open(fid, '/stage_movement2');
            if H5L.exists(gid,'stage_vec2','H5P_DEFAULT')
                H5L.delete(gid,'stage_vec2','H5P_DEFAULT');
            end
            
            if H5L.exists(gid,'is_stage_move2','H5P_DEFAULT')
                H5L.delete(gid,'is_stage_move2','H5P_DEFAULT');
            end
            
            if H5L.exists(gid,'frame_diff2','H5P_DEFAULT')
                H5L.delete(gid,'frame_diff2','H5P_DEFAULT');
            end
            H5L.delete(gid,'/stage_movement2','H5P_DEFAULT');
        end
        H5F.close(fid);
        
        
        %% save stage vector
        
        h5create(skeletons_file, '/stage_movement2/stage_vec2', size(stage_vec_d), 'Datatype', 'double', ...
            'Chunksize', size(stage_vec_d), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
        h5write(skeletons_file, '/stage_movement2/stage_vec2', stage_vec_d);
        
        h5create(skeletons_file, '/stage_movement2/is_stage_move2', size(is_stage_move_d), 'Datatype', 'int8', ...
            'Chunksize', size(is_stage_move_d), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
        h5write(skeletons_file, '/stage_movement2/is_stage_move2', is_stage_move_d);
        
        h5create(skeletons_file, '/stage_movement2/frame_diffs2', size(frame_diffs_d), 'Datatype', 'double', ...
            'Chunksize', size(frame_diffs_d), 'Deflate', 5, 'Fletcher32', true, 'Shuffle', true)
        h5write(skeletons_file, '/stage_movement2/frame_diffs2', frame_diffs_d);
        
        h5writeatt(skeletons_file, '/stage_movement2', 'fps', fps)
        h5writeatt(skeletons_file, '/stage_movement2', 'delay_frames', delay_frames)
        h5writeatt(skeletons_file, '/stage_movement2',  'pixel_per_micron_scale',  pixelPerMicronScale)
        h5writeatt(skeletons_file, '/stage_movement2',  'rotation_matrix',  rotation_matrix)
        
    end
end

%masked_image_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/MaskedVideos/goa-1 (sa734)I on food L_2010_03_04__10_44_32___8___6.hdf5';
%skeletons_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/Results/goa-1 (sa734)I on food L_2010_03_04__10_44_32___8___6_skeletons.hdf5';

