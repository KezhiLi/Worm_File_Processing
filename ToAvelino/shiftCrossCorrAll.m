%function alignStageMotion(masked_image_file,skeletons_file, is_swimming)
main_dir = '/Users/ajaver/Desktop/Videos/single_worm/agar_2/MaskedVideos/';
results_dir = strrep(main_dir, 'MaskedVideos', 'Results');


for file = dir(main_dir)'
    if ~isempty(regexp(file.name, '\w*.hdf5', 'ONCE'))
        disp(file.name)
        masked_image_file = fullfile(main_dir, file.name);
        skeletons_file = fullfile(results_dir, strrep(file.name, '.hdf5', '_skeletons.hdf5'));
        stage_file = fullfile(results_dir, strrep(file.name, '.hdf5', '_stagemov.mat'));
        
        video_timestamp_ind = h5read(skeletons_file, '/timestamp/raw');
        video_timestamp_time = h5read(skeletons_file, '/timestamp/time');
%%
        clear xShift yShift abs_diff_fra 
        
        [xShift, yShift, abs_diff_fra] = shiftCrossCorrelation(masked_image_file);
        
        
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

        save(stage_file, 'xShift', 'yShift', 'video_timestamp_ind', 'video_timestamp_time', 'stage_x', 'stage_y', 'stage_time')
    end
end


