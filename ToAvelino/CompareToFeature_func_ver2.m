%skeletons_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/Results/goa-1 (sa734)I on food L_2010_03_04__10_44_32___8___6_skeletons.hdf5';
%features_mat = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/Features/goa-1 (sa734)I on food L_2010_03_04__10_44_32___8___6_features.mat';
%masked_image_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/MaskedVideos/goa-1 (sa734)I on food L_2010_03_04__10_44_32___8___6.hdf5';

tt_step = 15;

%% segworm data
load(features_mat);
is_stage_move_seg = info.video.annotations.frames == 3;

ske_segworm_x = worm.posture.skeleton.x;
ske_segworm_y = worm.posture.skeleton.y;

%% read skeleton file and results from the stage alignment
%x_pixel_per_microns = 1/(h5readatt(skeletons_file, '/trajectories_data', 'pixels2microns_x'));
%y_pixel_per_microns = 1/(h5readatt(skeletons_file, '/trajectories_data', 'pixels2microns_y'));
x_pixel_per_microns = 1/(h5readatt(masked_image_file, '/mask', 'pixels2microns_x'));
y_pixel_per_microns = 1/(h5readatt(masked_image_file, '/mask', 'pixels2microns_y'));


stage_vec = h5read(skeletons_file, '/stage_vec');
skeleton_hdf5 = h5read(skeletons_file,'/skeleton');

n_segments = size(skeleton_hdf5,2);
frame_total = size(skeleton_hdf5,3);

video_timestamp_ind = int64(h5read(skeletons_file, '/timestamp/raw'));
video_timestamp_ind = video_timestamp_ind + 1; % correct for python indexing
tot_timestamp = max(video_timestamp_ind);
tot_timestamp = max(size(ske_segworm_x,2), tot_timestamp);

if tot_timestamp > size(ske_segworm_x,2)
    dd = tot_timestamp-size(ske_segworm_x,2);
    ske_segworm_x = padarray(ske_segworm_x, double([0, dd]), nan, 'post');
    ske_segworm_y = padarray(ske_segworm_y, double([0, dd]), nan, 'post');
end


% calculate skeleton coordinates after considering the stage moving
x_ske_cum = squeeze(skeleton_hdf5(1,:,:)) - ones(n_segments,1)*stage_vec(1,:);
y_ske_cum = squeeze(skeleton_hdf5(2,:,:)) - ones(n_segments,1)*stage_vec(2,:);


x_ske_aligned = nan(n_segments, tot_timestamp);
y_ske_aligned = nan(n_segments, tot_timestamp);
x_ske_aligned(:, video_timestamp_ind) = x_ske_cum*(x_pixel_per_microns);
y_ske_aligned(:, video_timestamp_ind) = y_ske_cum*(y_pixel_per_microns);


is_stage_move = isnan(stage_vec(1,:));
is_bad_skel = isnan(squeeze(skeleton_hdf5(1,1,:)));



x_ske_aligned = -x_ske_aligned;
dX = nanmedian(nanmedian(x_ske_aligned-ske_segworm_x));
x_ske_aligned = x_ske_aligned - dX;

y_ske_aligned = -y_ske_aligned;
dY = nanmedian(nanmedian(y_ske_aligned-ske_segworm_y));
y_ske_aligned = y_ske_aligned - dY;


figure, hold on
plot(x_ske_aligned(:, 1:tt_step:end)' , y_ske_aligned(:, 1:tt_step:end)', '.r')
plot(ske_segworm_x(:,1:tt_step:end)', ske_segworm_y(:, 1:tt_step:end)', '.b')
axis equal

%%
%{


%%
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

%%
fig66 = figure(66);
hold on

h1 = hggroup();
for tt_1 = 1:tt_step:size(worm_skeleton_x,2);
    if ~(isnan(sum(worm_skeleton_x(:,tt_1))) || isnan(sum(worm_skeleton_y(:,tt_1))))
        h = plot(worm_skeleton_x(:,tt_1), worm_skeleton_y(:,tt_1),'r' );
        set(h, 'parent', h1);
    end
end

h2 = hggroup();
plot(x_ske_cum_left2(1,:),y_ske_cum_left2(1,:),'b' ); axis equal; 
for tt_1 = 1:tt_step:size(x_ske_cum_left2,1);
    h = plot(x_ske_cum_left2(tt_1,:),y_ske_cum_left2(tt_1,:) ,'b');
    set(h, 'parent', h2);
end
legend([h1, h2], 'Segworm', 'New')


frame_total_ske = sum(~isnan(x_ske_cum_left(:,1)));
% show frames shown
frame_show1 = size(x_ske_cum_left2,1);

%%

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
peak_ind_in_total(video_timestamp_ind(cancel_fra_ind2)+1) = 1.15;
%hold on, plot(peak_ind_in_total,'r')
hold on, plot(peak_ind_in_total,'b')
hold off,

% stage_ind = zeros(frame_total,1);
% stage_ind(cancel_fra_ind2+1) = 1;
% if stage_ind(2)==1
%     stage_ind(1) = 1;
% end

% stage_vec_save = [trajectories_folder,name,'_stage_vec.mat'];
% save(stage_vec_save,'stage_ind');


%
%
% frame_all = zeros(1, size( worm.posture.skeleton.x,2));
% frame_all1 = frame_all;
% frame_all1(diff_mask_central(cancel_fra_ind,8)) = 1;
% frame_stage = frame_all;
% frame_stage(info.video.annotations.frames==3) = 1.3;
%
% fig71 = figure(71), plot(frame_all1,'g');
% hold on, plot(frame_stage,'r')
%}