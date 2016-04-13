# Worm_File_Processing 

**Objective:**These files are for the alignement of worm vedios. Because the camera stages moves when the worm moves out of the central of scope, we need to build a correspondence between the stage motion recorded and the real stage motion time in the video.

The Worm Alignement consists of 4 steps:

1. **Read Stage Motion Data:** There are two ways to read the necessary data. 

*Note: the tracker uses numpy arrays with the C ordering (the last dimension is the fast changing).*


  - `/mask` compressed array with the masked image.

