
%function alignStageMotion(masked_image_file,skeletons_file, is_swimming)
main_dir = 'Z:\single_worm\agar_2\MaskedVideos\';
results_dir = strrep(main_dir, 'MaskedVideos', 'Results');
feat_dir = strrep(main_dir, 'MaskedVideos', 'Features');
is_swimming = false;

%for file = dir(main_dir)'
    if ~isempty(regexp(file.name, '\w*.hdf5', 'ONCE'))
        masked_image_file = fullfile(main_dir, file.name);
        skeletons_file = fullfile(results_dir, strrep(file.name, '.hdf5', '_skeletons.hdf5'));
        features_mat = fullfile(feat_dir, strrep(file.name, '.hdf5', '_features.mat'));
        %%
        fid = H5F.open(skeletons_file,'H5F_ACC_RDWR','H5P_DEFAULT');
        is_finished = H5L.exists(fid,'stage_vec','H5P_DEFAULT');
        H5F.close(fid);
        %%
        %disp(file.name)
        if exist(features_mat, 'file') % && is_finished
            
            disp(features_mat)
            CompareToFeature_func_ver2
            %break
        end
        %

    end
%end


