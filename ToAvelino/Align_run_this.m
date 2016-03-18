clear
clc

%% function input parameters

%%Input: 
name = 'egl-23 (n601)IV on food R_2010_07_20__12_24_56___1___7';
trajectories_folder = 'Z:\Results_old\nas207-3\Data\from pc207-13\Laura\20-07-10\4\';
trajectories_file = [trajectories_folder,name,'_trajectories.hdf5'];
hdf5_path =['Z:\MaskedVideos_old\nas207-3\Data\from pc207-13\Laura\20-07-10\4\',name,'.hdf5'];
ske_file =[trajectories_folder,name,'_skeletons.hdf5'];
xml_file = ['Z:\thecus\nas207-3\Data\from pc207-13\Laura\20-07-10\4\.data\',name,'.info.xml'];
csv_file =['Z:\thecus\nas207-3\Data\from pc207-13\Laura\20-07-10\4\.data\',name,'.log.csv'];

features_mat ='N:\Andre\results-12-05-10\Laura Grundy\egl-23\n601\MT1231\on_food\XX\30m_wait\L\tracker_1\2010-07-20___12_24_56\egl-23 (n601)IV on food R_2010_07_20__12_24_56___1___7_features.mat';

%% Outout

[stage_vec,pixel_to_micro] = Align_TimeStamp_Func(name,trajectories_folder,trajectories_file,hdf5_path,ske_file,xml_file,csv_file,features_mat);