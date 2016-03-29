stage_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_1/Results/431 JU298 on food L_2011_03_17__12_02_58___2___3_stagemov.mat';
load(stage_file)
%%

stage_x_n = stage_x-stage_x(1);
stage_y_n = stage_y-stage_y(1);

%%
%it seems that the stage and the worm coord are shifted...
img_shift_x = [0; cumsum(xShift)];
img_shift_y = [0; cumsum(yShift)];

%%
figure, hold on
plot(stage_time, stage_x_n, 'color', [0.8500    0.3250    0.0980])
plot(video_timestamp_time, img_shift_x, 'r')

figure, hold on
plot(video_timestamp_time, img_shift_y, 'b')
plot(stage_time, stage_y_n)
%%
stage_dx = diff(stage_x);
stage_dy = diff(stage_y);
stage_dt = diff(stage_time);



movie_period = median(diff(video_timestamp_time));
%check there are not elements repeated
assert(numel(unique(stage_time)) == numel(stage_time));

stage_min_delay = min(stage_dt);
stage_window = round(stage_min_delay/movie_period);
%%
di = stage_window;
dx = img_shift_x(di+1:end) - img_shift_x(1:end-di);


figure, hold on
plot(video_timestamp_time(1:end-di), dx)
plot(stage_time(1:end-1), stage_dx, 'o')
%%
labels_x = zeros(size(xShift));
group_n = 0;
is_group = false;
for ii = 1:numel(xShift);
    if xShift(ii) ~= 0
        if ~is_group
            group_n = group_n + 1;
            is_group = true;
        end
        labels_x(ii) = group_n;
    else
        is_group = false; 
    end
end

tot_groups = numel(unique(labels_x)) - 1;
x_img_dist = zeros(1, tot_groups);

for ii = 1:numel(xShift);
    if labels_x(ii) > 0
        group_n = labels_x(ii);
        x_img_dist(group_n) = x_img_dist(group_n) + xShift(ii);
    end
end
%%
figure, hold on
plot(x_img_dist)
plot(stage_dx)


%%




%%
%{
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
%}
