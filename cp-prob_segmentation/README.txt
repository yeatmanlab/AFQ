SEGMENTATION OF THE CEREBELLAR PEDUNCLES IN PROBABILISTIC TRACTOGRAPHY 
-----------------------------------------------------------------------

This directory contains functions, scripts, and ROIs to identify the human cerebellum-cerebrum white matter connections, namely the inferior (ICP), middle (MCP) and superior (SCP) cerebellar peduncles.

Folders: ROIs - Contains the region of interests (ROIs) in nifti format used to segment the cerebellar peduncles. Note: the folder can be saved anywhere since the code will ask you to choose the right directory

Functions: AFQ_SegmentCerebellum - Main function that takes an existing afq structure (default 20 fiber groups) and segments and adds the cerebellar peduncles, i.e. left/right ICP, left/right SCP and MCP.

Additional scripts - Modified version of existing AFQ functions called by AFQ_SegmentCerebellum. Functions need to be saved to in the AFQ/function directory and added to your MATLAB path. Note: Function names have been changed to begin with "cp_" so that they don't mask any existing AFQ functions.

If you use this code for your own study, please cite the following article as a reference:

S. Jossinger, M. Yablonski, O. Amir, M. Ben-Shachar; The contributions of the cerebellar peduncles and the frontal aslant tract in mediating speech fluency. Neurobiology of Language 2023; doi: https://doi.org/10.1162/nol_a_00098