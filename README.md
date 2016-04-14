# Worm_File_Processing 

**Objective:**These files are for the alignement of worm vedios. Because the camera stages moves when the worm moves out of the central of scope, we need to build a correspondence between the stage motion recorded and the real stage motion time in the video.

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

* Note: The following scripts are based the result of 1), given (`masked_image_file`, `skeletons_file`).

2. **perform improved alignment algorithm**
This improved alignment algorithm is based on segworm, improved by Avelino Javer and Kezhi Li. The idea of segworm's alignment is based on the variance of the pixel difference between succesive frames (`frame_diffs`).  



