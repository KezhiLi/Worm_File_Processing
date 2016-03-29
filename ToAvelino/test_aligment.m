%masked_image_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/MaskedVideos/goa-1 (sa734)I on food R_2009_07_08__12_09__6.hdf5';
%skeletons_file = '/Users/ajaver/Desktop/Videos/single_worm/agar_goa/Results/goa-1 (sa734)I on food R_2009_07_08__12_09__6_skeletons.hdf5';
%is_swimming = false;

alignStageMotion(masked_image_file,skeletons_file, is_swimming);

