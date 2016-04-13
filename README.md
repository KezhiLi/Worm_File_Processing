# Worm_File_Processing 

**Objective:**These files are for the alignement of worm vedios. Because the camera stages moves when the worm moves out of the central of scope, we need to build a correspondence between the stage motion recorded and the real stage motion time in the video.

The Worm Alignement consists of 4 steps:

1. **Read Stage Motion Data:** There are two ways to read the necessary data: from 1) .hdf5 file or 2) .xml, .csv and other files
1) Two .hdf5 files (`masked_image_file`, `skeletons_file`) include all necessary data to proceed the alignment, including
  - `/mask` compressed array with the masked image.
  - `/xml_info`


*Note: the tracker uses numpy arrays with the C ordering (the last dimension is the fast changing).*



