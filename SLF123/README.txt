SEGMENTATION OF THE THREE SLF BRANCHES
-------------------------------------------------------------------------------------

This directory contains functions and regions of interest (ROIs) to segment the human
superior longitudinal fasciculus (SLF) into three branches: SLF-I (dorsal branch), 
SLF-II (middle branch), and SLF-III (ventral branch), automatically using AFQ.
New ROIs were defined on the MNI T2 template in accordance with defenitions in 
Thiebaut de Schotten et al. (2011) Nat. Neurosci. paper.

ROIs folder: contains ROIs in nifti format.
Note: the folder can be saved anywhere since the code will ask you to choose the 
directory where they are stored.

Each SLF branch is defined by three ROIs: frontal and parietal 'AND' ROIs and a 'NOT' 
temporal ROI to exclude frontotemporal fibers of the arcuate fasciculus. Frontal ROIs 
are located on the coronal plane of the AC and defined distinctly for each of the 
three branches (SLF-I - 'SFgL'/'R', SLF-II - 'MFgL'/'R', SLF-III - 'PrgL'/'R'). The 
parietal ROI ('PaL'/'R') is shared for the three branches, located on the coronal 
plane of the PC. The temporal 'NOT' ROI is also shared, and is identical to the AFQ 
ROI2 used to segment the arcuate fasciculus ('SLFt_roi2_L'/'R').

Functions: 
AFQ_SegmentSLF123 - main function that will add the left and right SLF-I,-II, and -III
to an existing afq structure (default 20 mori fiber groups). The tracts will be added 
as fiber groups number 21-26.

SLF123_AFQ_AddNewFiberGroup3ROI - a modified version of AFQ_AddNewFiberGroup. The 
main function will call it to segment the bilateral three SLF branches and compute 
their tract properties.

For details on preprocessing params see Methods in the attached manuscript 'White 
matter associations with spelling performance' (Sagi et al, under review) 

Romi Sagi, August 2023