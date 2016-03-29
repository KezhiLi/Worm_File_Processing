%%
masked_image_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_1/MaskedVideos/431 JU298 on food L_2011_03_17__12_02_58___2___3.hdf5';
[masked_dir, fname, ext] = fileparts(masked_image_file);
results_dir = strrep(masked_dir, 'MaskedVideos', 'Results');
feat_dir = strrep(masked_dir, 'MaskedVideos', 'Features');

skeletons_file = fullfile(results_dir, [fname, '_skeletons.hdf5']);
features_mat = fullfile(feat_dir, [fname, '_features.mat']);
        
%% read pixels per microns attribute. This info was extracted from the info.xml file
x_pixel_per_microns = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_x');
y_pixel_per_microns = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_y');
%% read stage data from the mask hdf5 file. This information was extracted
%from the .log.csv file 
stage_data = h5read(masked_image_file, '/stage_data');
%correct for duplicated data keeping the last instance of a given time
[stage_time, ind] = unique(stage_data.stage_time, 'last');
%stage_time = stage_time*60; %Convert the time in seconds. It was originally in minutes.
stage_x = stage_data.stage_x(ind)/x_pixel_per_microns;
stage_y = stage_data.stage_y(ind)/y_pixel_per_microns;

stage_x_n = stage_x-stage_x(1);
stage_y_n = stage_y-stage_y(1);

%% the centroid position and real time frame it's better stored in the skeletons_file

video_timestamp_ind = h5read(skeletons_file, '/timestamp/raw');
video_timestamp_time = h5read(skeletons_file, '/timestamp/time');
trajectories_data = h5read(skeletons_file, '/trajectories_data');

%%

skeletons = h5read(skeletons_file, '/skeleton');

x_diff = squeeze(diff(skeletons(1,:,:), 1, 3));
y_diff = squeeze(diff(skeletons(2,:,:), 1, 3));
%% get a shift form the skeletons
[~, dx_imin] = min(abs(x_diff));
[~, dy_imin] = min(abs(y_diff));

good_ind_x = sub2ind(size(x_diff), dx_imin, 1:size(x_diff,2));
dx_min = x_diff(good_ind_x);

good_ind_y = sub2ind(size(y_diff), dx_imin, 1:size(y_diff,2));
dy_min = y_diff(good_ind_y);


good = ~(isnan(dy_min) | isnan(dx_min));
shift_skel_y = cumsum(dy_min(good));
shift_skel_x = cumsum(dx_min(good));
shift_skel_t = trajectories_data.timestamp_time(good);

figure, hold on
plot(stage_time, stage_x_n)
plot(shift_skel_t, shift_skel_x, 'm')
%% calculate shift from cross correlation between frames, and get the absolute difference between images
[xShift, yShift, abs_diff_fra] = shiftCrossCorrelation(masked_image_file);
%%
%it seems that the stage and the worm coord are shifted...
dd = yShift;
yShift = xShift;
xShift = dd;

img_shift_x = [0; cumsum(xShift)];
img_shift_y = [0; cumsum(yShift)];

%%
figure, hold on
plot(stage_time, stage_x_n, 'color', [0.8500    0.3250    0.0980])
plot(video_timestamp_time, img_shift_x, 'r')

figure, hold on
plot(video_timestamp_time, img_shift_y, 'b')
plot(stage_time, stage_y_n)

%% resample to be able to get a higher resolution in the shift space
deltaT = median(diff(video_timestamp_time));
max_time = max(max(video_timestamp_time), max(stage_time));
tot_resampled = round(max_time/deltaT);

extra_shift = 15;

img_shift_x_pad = padarray(img_shift_x, extra_shift, 'replicate');
img_shift_y_pad = padarray(img_shift_y, extra_shift, 'replicate');

time_dilations = 0.75:0.01:1.25;

for kk = 1%:numel(time_dilations)
    current_shift = extra_shift;
    t_dil = 1;%time_dilations(kk);
    
    stage_ind = round(stage_time*t_dil/deltaT) + current_shift;
    
    x_stage_impulse = zeros(1, tot_resampled + extra_shift*2);
    x_stage_impulse(stage_ind) = stage_x_n;
    
    %fill with the next value from: http://stackoverflow.com/questions/7348615/matlab-fill-previous-value-if-missing-value-or-zero
    x_stage_impulse(stage_ind(2:end)) = x_stage_impulse(stage_ind(2:end)) - x_stage_impulse(stage_ind(1:end-1));
    x_stage_full = cumsum(x_stage_impulse);
    
    
    figure, hold on
    plot(x_stage_full)
    plot(img_shift_x_pad)
end
%%



%%
figure, hold on
plot(trajectories_data.timestamp_time, trajectories_data.cnt_coord_x) 
plot(trajectories_data.timestamp_time, trajectories_data.cnt_coord_y)
%%
figure, hold on
plot(trajectories_data.timestamp_time, trajectories_data.coord_x) 
plot(trajectories_data.timestamp_time, trajectories_data.coord_y)
%%
plate_worms = h5read(trajectories_file, '/plate_worms');
good = plate_worms.worm_index_joined==1;
figure, hold on
plot(plate_worms.coord_x(good)) 
plot(plate_worms.coord_y(good))


%%
%{
delta_t = 0;
yy = trajectories_data.cnt_coord_y;
dy = yy(1:end-delta_t) - yy(delta_t+1:end);
tt = trajectories_data.timestamp_time(1:end-delta_t);
figure, hold on
plot(tt, cumsum(dy))
plot(stage_time, stage_y-stage_y(1))
%}