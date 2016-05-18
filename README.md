# Worm_File_Processing 

**Objective:**This package is for the alignement of worm vedios. Because the camera stages moves when the worm moves out of the central of the scope, we need to build a correspondence between the stage motion recorded and the real stage motion time in the video. Modify the file `ReadFiles.m` accordingly to read file names that are going to process.  In it, the function ` alignStageMotionFun.m` is run to do the Worm Alignment. The inputs are `.hdf5` and `_skeletons.hdf5`. The alignment results are saved in `_skeletons.hdf5` with fields `/stage_movement`.

The Worm Alignement consists of 4 steps:

1. **Read Stage Motion Data:** There are two ways to read the necessary data: from 1) .hdf5 file or 2) .xml, .csv and other files 

  1) Two .hdf5 files (`masked_image_file`, `skeletons_file`) include all necessary data to proceed the alignment, including
  - `/mask` compressed array with the masked image, and pixel per microns.
  - `/xml_info` delay_frames
  - `/timestamp/time` the time stamp, which can lead to fps
  - `/trajectories_data` calculate the central mask shift
  
  2) `.xml` reads pixel per microns; `.csv` reads the stage motion postions and time; `.hdf5` reads `/mask`; `_trajectories.hdf5` reads
  - `/timestamp`  time stamp index
  - `/timestamp_time`  time stamp time

  Note: The following scripts are based the result of 1), given (`masked_image_file`, `skeletons_file`).

2. **perform basic alignment algorithm** This basic alignment algorithm is based on segworm alignment, improved by Avelino Javer and Kezhi Li. The idea of segworm's alignment is based on the variance of the pixel difference between succesive frames (`frame_diffs` is calculated by `getFrameDiffVar_fast.m`, a 10X fast version of original function) and identify peaks one by one (this step is run in `findStageMovement_gs.m` when `wind_weights`=0).  For most cases, it will obtain the exact same result as segworm alignment.

 ![frame_diffs](https://github.com/KezhiLi/Worm_File_Processing/blob/master/frame_diffs_github1.png?raw=true)
 
   Note: 
   
   1) Avelino modified the segmentation algorithm (set threshold adaptively and smoothly), refined the skeleton and read `.hdf5` files instead of `.avi` or `.mpeg` files comparing to segworm. Kezhi continues to update the algorithm in the aspects of
    - after recognizing 1 peak, set it to 0 to prevent the miss recognition of next peak
    - switch the start/end frame index when the later is smaller than the fronter
    - adjust `otsuThr` if more than 1 peak is found in the frames interval in consideration
   
   2) before running `findStageMovement_gs.m`, in `alignStageMotionFun.m` some preparation is done, such as deleting `stage_vec`, `is_stage_move`, `frame_diffs` if there exist (because they the results will be saved).   

  If the alignment process goes smoothly, the function `findStageMovement_gs.m` outputs `[is_stage_move, movesI, stage_locations]` as results. Specifically
  - `is_stage_move` a binary vector that shows 1 when a frame belongs to stage motion and 0 otherwise
  - `movesI`: an `N*2` vector that indicates the start/end frame index of each stage motion(peak), `N` is the number of motions 
  - `stage_location` an `N*2` vector to show the x-pixel/y-pixel location after compensating the stage motion
  
  Finally, these results will be converted to `is_stage_move_d` and `stage_vec_d` and saved in `_skeletons.hdf5` in field `/stage_movement`.

3. **perform advanced alignment algorithm** This advanced alignment algorithm is based on segworm alignment, and modified by generate a Gaussian window over the frames interval in consideration. It is performed in `alignStageMotionFun.m` when step 2 fails. The new stage motions are esimated by `findStageMovement_gs.m` when `wind_weights`>0, and uses a new `frame_diffs`.

  The new `frame_diffs` is masked by a logical operation of tow criterions.
  - `diff_mask_central_abs` the absolute distance of mask central shifts
  - `[xShift, yShift]` the cross-correlation of the pixels of mask
  After thresholding, imdilating, filtering, fixing (including filling and trimming), we obtain the new `frame_diffs` as one input to 
  `findStageMovement_gs.m`. 

  In `findStageMovement_gs.m`, a Gaussin window is generated with parameter `wind_weights`. It implies that the peaks that near the predicted time is strengthened, while peaks far from the predicted frame time will be suppressed. This operation is to reduce the affect of wrong peaks that happens before the right time (in that case the wrong peak could be recognized by mistake). 
  
   Note:

   1) `wind_weights` is a parameter to determine the arch of Gaussian window. When `wind_weights`=0, the algorithm returns to algorithm similar to conventional segworm alignment; When `wind_weights` =1, the algorithm turns to one with high ceiling Gaussian window, which pays more attention to the central. The stratigy in step 3 is increasing `wind_weights` incrementally to expect a outcomes that has minimal modication of the alignment function.

4. **check the result and save** The alighment results are saved in `_skeletons.hdf5` in field `/stage_movement` with different `exit_flag`. In detail, the following vectors are saved
  - `frame_diffs`: the differences between frames (can be old or updated)
  - `stage_vec` vector of x,y locations after alignment adjustment to compensate the stage motion
  - `is_stage_move` a binary vector to tell if one frame belongs stage motion, the same as explained in step 2
  - `has_finished` it is `exit_flag`, means
      - `1` find stage motion smoothly
      - `80` the timestamp is corrupted or does not exist
      - `81` number of timestamps does not match the number read movie frames
      - `82` cannot find stage motion smoothly in the basic algorithm, then it can be changed to other numbers:
        - `2` find stage motion smoothly in the advanced algorithm
        - `71` hdf5 and skeleton has different numbers of frames
        - `72` still cannot find stage motion in the advanced algorithm
  - `fps` frame rate
  - `delay_frames` most number of frames that a delay can happen
  - `pixel_per_micron_scale` pixel per micron that can be used for converting pixel value to microns value and vise versa.
  - `rotation_matrix` compensate the rotation error in calibration

**Interface that adjust the alignment manually:** we also generate a user-interface to users to adjust the alignment manually when the automatic alignment step fails. The main interface looks like the one below, and the key panel is 'Manually moidfy frameDiffs'.

 ![frame_diffs](https://github.com/KezhiLi/Worm_File_Processing/blob/master/Manual_interface.png)
