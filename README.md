# Worm_File_Processing 

**Objective:**These files are for the alignement of worm vedios. Because the camera stages moves when the worm moves out of the central of scope, we need to build a correspondence between the stage motion recorded and the real stage motion time in the video. Modify the file `ReadFiles.m` accordingly to read file names that are going to process.  In it, the function ` alignStageMotionFun.m` is run to do the Worm Alignment. The inputs are `.hdf5` and `_skeletons.hdf5`, and the alignment results are saved in `_skeletons.hdf5` with fields `/stage_movement`.

The Worm Alignement consists of 4 steps:

1. **Read Stage Motion Data:** There are two ways to read the necessary data: from 1) .hdf5 file or 2) .xml, .csv and other files 

  1) Two .hdf5 files (`masked_image_file`, `skeletons_file`) include all necessary data to proceed the alignment, including
  - `/mask` compressed array with the masked image, and pixel per microns.
  - `/xml_info` delay_frames
  - `/timestamp/time` the time stamp, which can lead to fps
  - `/trajectories_data` calculate the central mask shift
  
  2) `.xml` reads pixel per microns, `.csv` reads the stage motion postions and time, `.hdf5` reads `/mask`, `_trajectories.hdf5` reads
  - `/timestamp`  time stamp index
  - `/timestamp_time`  time stamp time

  Note: The following scripts are based the result of 1), given (`masked_image_file`, `skeletons_file`).

2. **perform basic alignment algorithm** This basic alignment algorithm is based on segworm alignment, improved by Avelino Javer and Kezhi Li. The idea of segworm's alignment is based on the variance of the pixel difference between succesive frames (`frame_diffs`) and identify peaks one by one (this is run in `findStageMovement_gs.m` when `wind_weights`=0).  

 ![frame_diffs](https://github.com/KezhiLi/Worm_File_Processing/blob/master/ToAvelino/stage_correction_segworm/frame_diffs_github1.png?raw=true)
 
   Note: 
   1) Avelino modified the segmentation algorithm (set threshold adaptively and smoothly), refined the skeleton and read `.hdf5` files instead of `.avi` or `.mpeg` files comparing to segworm. Kezhi continues to update the algorithm in the aspects of
      - after recognizing 1 peak, set it to 0 to prevent the miss recognition of next peak
      - switch the start/end frame index when the later is smaller than the fronter
      - adjust `otsuThr` if more than 1 peak is found in the frames interval in consideration
   2) before running `findStageMovement_gs.m`, in `alignStageMotionFun.m` some preparation is done, such as delete `stage_vec`, `is_stage_move`, `frame_diffs` if there exist (because they the results will be saved).   




