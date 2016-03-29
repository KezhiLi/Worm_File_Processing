%function alignStageMotion(masked_image_file,skeletons_file, is_swimming)

main_dir = '/Users/ajaver/Desktop/Videos/single_worm/agar_2/MaskedVideos/';
results_dir = strrep(main_dir, 'MaskedVideos', 'Results');
feat_dir = strrep(main_dir, 'MaskedVideos', 'Features');

is_swimming = false;

for file = dir(main_dir)'
    if ~isempty(regexp(file.name, '\w*.hdf5', 'ONCE'))
        disp(file.name)
        masked_image_file = fullfile(main_dir, file.name);
        skeletons_file = fullfile(results_dir, strrep(file.name, '.hdf5', '_skeletons.hdf5'));
        features_mat = fullfile(feat_dir, strrep(file.name, '.hdf5', '_features.mat'));
        
        fid = H5F.open(skeletons_file,'H5F_ACC_RDWR','H5P_DEFAULT');
        is_finished = H5L.exists(fid,'stage_vec','H5P_DEFAULT');
        H5F.close(fid);
        
        if ~is_finished && exist(features_mat, 'file')
            [stage_vec, pixel_to_micro] = Align_TimeStamp_Func_ver2(masked_image_file,skeletons_file, is_swimming);
        end
        
    end
end

%masked_image_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/MaskedVideos/goa-1 (sa734)I on food L_2010_03_04__10_44_32___8___6.hdf5';
%skeletons_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/Results/goa-1 (sa734)I on food L_2010_03_04__10_44_32___8___6_skeletons.hdf5';

