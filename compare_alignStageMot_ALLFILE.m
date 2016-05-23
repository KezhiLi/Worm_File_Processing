%function: see the alignment peaks and skeletons of the recovered hdf5
%files

clear
clc

iif = 18;

% change '/' to '\' due to the difference between python and matlab
failed_files_all = strrep(fileread('good_files.txt'),'/','\');

% replace folder
gap_sym = '\Volumes\behavgenom_archive$';

ini_loc = strfind(failed_files_all,gap_sym);
%ini_loc = regexp(failed_files_all,gap_sym);

file_name = {};

% restore file names to independent cell
for ii = 1:numel(ini_loc)-1
    file_name = [file_name;failed_files_all(ini_loc(ii):ini_loc(ii+1)-2)];
end
file_name = [file_name;failed_files_all(ini_loc(numel(ini_loc)):end)];

% set current file and result hdf5 file
cur_file = strtrim(file_name{iif});

cur_file = strrep(cur_file, '\Volumes\behavgenom_archive$\', 'Z:\');

masked_image_file = cur_file;

skeletons_file = strrep(strrep(cur_file, '.hdf5', '_skeletons.hdf5'),'MaskedVideos','Results_old');
% features_mat = fullfile(strrep(cur_file, '.hdf5', '_features.mat'));

time_ske = h5read(skeletons_file,'/timestamp/time');
stage_data = h5read(masked_image_file, '/stage_log');
mediaTimes = stage_data.stage_time';%*60;
ii =1;
while ii == 1;
    if mediaTimes(ii:ii+1) == [0,0];
        mediaTimes = mediaTimes(2:end);
    else
        break;
    end
end

record_stage = zeros(size(time_ske));
for ii = 2:length(mediaTimes)
    [minn, min_loc] =min(abs(mediaTimes(ii)-time_ske));
    record_stage(min_loc) = 1;
end




is_swimming = false;


%stage_vec1 = h5read(skeletons_file, '/stage_movement/stage_vec');
stage_vec2 = h5read(skeletons_file, '/stage_movement/stage_vec');

%is_stage_move1 = h5read(skeletons_file, '/stage_movement/is_stage_move');
is_stage_move2 = h5read(skeletons_file, '/stage_movement/is_stage_move');

rotation_matrix = h5readatt(skeletons_file, '/stage_movement' ,'rotation_matrix');
pixelPerMicronScale = h5readatt(skeletons_file, '/stage_movement','pixel_per_micron_scale');
frame_diffs_d = h5read(skeletons_file, '/stage_movement/frame_diffs');

%stage_vec_abs = sqrt(stage_vec1(1,:).^2+stage_vec1(2,:).^2);
%stage_vec2_abs = sqrt(stage_vec2(1,:).^2+stage_vec2(2,:).^2);


figure,plot([0,frame_diffs_d/max(frame_diffs_d)*4]); 
hold on, plot( is_stage_move2,'c');
plot(record_stage*0.5,'r');

ind_peak = find(record_stage>0.5);
cur_p = 1;
for ii =1:length(is_stage_move2)-1
    if is_stage_move2(ii:ii+1) == [0,1];
        plot([ii+1,ind_peak(cur_p)],[1.5,0.5], 'm');
        cur_p = cur_p +1;
    end
end
hold off

%vec_diff1 =stage_vec1(1,:)-stage_vec2(1,:);
%vec_diff2 =stage_vec1(2,:)-stage_vec2(2,:);

%disp([sum(abs(vec_diff1(~isnan(vec_diff1)))),sum(abs(vec_diff2(~isnan(vec_diff2))))])
%disp([sum(is_stage_move1),sum(is_stage_move2)])

%figure,plot(is_stage_move1 - is_stage_move2)

%%
%
is_stage_move_d = is_stage_move2;
stage_vec_d = stage_vec2;

% load(features_mat)
% seg_motion = info.video.annotations.frames==2;
% worm_ske_x = worm.posture.skeleton.x(1,:);
% worm_x_ind = find(~isnan(worm_ske_x)==1);
% worm_ske_x_real = worm_ske_x(worm_x_ind);
% figure,plot(worm.posture.skeleton.x(:, 1:15:end),worm.posture.skeleton.y(:, 1:15:end))
% axis equal
% %         if (all(seg_motion==is_stage_move_d))
% %             disp('Segworm and this code have the same frame aligment.')
% %         end


skeletons = h5read(skeletons_file, '/skeleton');

skeletons_mu = nan(size(skeletons));

for kk = 1:size(skeletons_mu,3)-15
    pixels = skeletons(:, :, kk)';
    origin = stage_vec_d(:,kk);
    % Rotate the pixels.
    pixels = (rotation_matrix * pixels')';
    
    % Convert the pixels coordinates to micron locations.
    microns(:,1) = origin(1) - pixels(:,1) * pixelPerMicronScale(1);
    microns(:,2) = origin(2) - pixels(:,2) * pixelPerMicronScale(2);
    %             microns(:,1) = origin(1) + pixels(:,1) * pixelPerMicronScale(1);
    %             microns(:,2) = origin(2) + pixels(:,2) * pixelPerMicronScale(2);
    skeletons_mu(:,:,kk) = microns';
end


skel_x = squeeze(skeletons_mu(1,:,:));
skel_y = squeeze(skeletons_mu(2,:,:));

%         figure, plot(worm.posture.skeleton.x(25,:));
%         hold on, plot(skel_x(25,1:400))
%
%         figure
%         plot(squeeze(skel_x(25,:)))

figure
plot(skel_x(:, 1:15:end), skel_y(:, 1:15:end))
axis equal

figure
plot(skel_x(24, :), skel_y(24, :))
axis equal
% disp([sum(~isnan(worm.posture.skeleton.x(1,:))),sum(~isnan(skel_x(1,:)))])