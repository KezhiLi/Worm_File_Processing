%function alignStageMotion(masked_image_file,skeletons_file, is_swimming)


main_dir = 'Z:\single_worm\agar_2\MaskedVideos';
results_dir = strrep(main_dir, 'MaskedVideos', 'Results');
feat_dir = strrep(main_dir, 'MaskedVideos', 'Features');

addpath(genpath([results_dir,'.']));

is_swimming = false;



files = dir(main_dir);
iif = 4; %1:numel(files)


file = files(iif);

fprintf('%i) %s\n', iif, file.name)
masked_image_file = fullfile(main_dir, file.name);
skeletons_file = fullfile(results_dir, strrep(file.name, '.hdf5', '_skeletons.hdf5'));
features_mat = fullfile(feat_dir, strrep(file.name, '.hdf5', '_features.mat'));

%stage_vec1 = h5read(skeletons_file, '/stage_movement/stage_vec');
stage_vec2 = h5read(skeletons_file, '/stage_movement2/stage_vec2');

%is_stage_move1 = h5read(skeletons_file, '/stage_movement/is_stage_move');
is_stage_move2 = h5read(skeletons_file, '/stage_movement2/is_stage_move2');

rotation_matrix = h5readatt(skeletons_file, '/stage_movement2' ,'rotation_matrix');
pixelPerMicronScale = h5readatt(skeletons_file, '/stage_movement2','pixel_per_micron_scale');
frame_diffs_d = h5read(skeletons_file, '/stage_movement2/frame_diffs2');

%stage_vec_abs = sqrt(stage_vec1(1,:).^2+stage_vec1(2,:).^2);
%stage_vec2_abs = sqrt(stage_vec2(1,:).^2+stage_vec2(2,:).^2);

figure,plot([0,frame_diffs_d/max(frame_diffs_d)*4]); hold on, plot( is_stage_move2);

%vec_diff1 =stage_vec1(1,:)-stage_vec2(1,:);
%vec_diff2 =stage_vec1(2,:)-stage_vec2(2,:);

%disp([sum(abs(vec_diff1(~isnan(vec_diff1)))),sum(abs(vec_diff2(~isnan(vec_diff2))))])
%disp([sum(is_stage_move1),sum(is_stage_move2)])

%figure,plot(is_stage_move1 - is_stage_move2)

%%
%
is_stage_move_d = is_stage_move2;
stage_vec_d = stage_vec2;

load(features_mat)
seg_motion = info.video.annotations.frames==2;
worm_ske_x = worm.posture.skeleton.x(1,:);
worm_x_ind = find(~isnan(worm_ske_x)==1);
worm_ske_x_real = worm_ske_x(worm_x_ind);
figure,plot(worm.posture.skeleton.x(:, 1:15:end),worm.posture.skeleton.y(:, 1:15:end))
axis equal
%         if (all(seg_motion==is_stage_move_d))
%             disp('Segworm and this code have the same frame aligment.')
%         end


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

disp([sum(~isnan(worm.posture.skeleton.x(1,:))),sum(~isnan(skel_x(1,:)))])