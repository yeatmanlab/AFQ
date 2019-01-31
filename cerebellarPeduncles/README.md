
CP (version 1.0)

============================================
SEGMENTATION OF THE CEREBELLAR PEDUNCLES
============================================

This directory contains functions, scripts, and ROIs to identify the human cerebellum-cerebrum white matter connections, namely the inferior (ICP), middle (MCP) and superior (SCP) cerebellar peduncles.

Folders:
ROIs - Contains the region of interests (ROIs) in nifti format used to segment the cerebellar peduncles. Note: the folder can be saved anywhere since the code will ask you to choose the right directory

Functions:
AFQ_SegmentCerebellum - Main function that takes an existing afq structure (default 20 fiber groups) and segments and adds the cerebellar peduncles, i.e. left/right ICP, left/right SCP and MCP.

ClipFibers - Saves a separate mat file containing the clipped cerebellar peduncles. Note: This function is currently called as a part of AFQ_SegmentCerebellum -- it very useful for visual inspection of individual fiber renderings to ensure that each peduncle conforms to anatomical norms, but it's not necessary. 

Additional scripts - Modified version of existing AFQ functions called by AFQ_SegmentCerebellum. Functions need to be saved to in the AFQ/function directory and added to your MATLAB path. Note: Function names have been changed to begin with "cp_" so that they don't mask any existing AFQ functions. 


If you use this code for your own study, please cite the following article as a reference:  

Bruckert L, Shpanskaya K, McKenna ES, Borchers LR, Yablonski M, Blecher T, Ben-Schachar M, Travis KE, Feldman HM. Age-Dependent White Matter Characteristics of the Cerebellar Peduncles from Infancy Through Adolescence. The Cerebellum. 2019;1–16. Available from: http://link.springer.com/10.1007/s12311-018-1003-9

