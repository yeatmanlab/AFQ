#mrQ manual#
mrQ is a software package designed to calculate MR parameters (T1 and PD) using spoiled gradient echo scans (SPGR, FLASH). mrQ allows the evaluation of macromolecule tissue fraction (MTV) and the volume of interaction proton (VIP) as well as the surface interaction ratio (SIR). 

The software and the tissue parameters are describe in the following article

>Mezer A, Yeatman JD, Stikov N, Kay K, Cho NJ, Dougherty R, Perry LM, Parvizi J, Hua L, Butts-Pauly K, Wandell BA. Measuring within the voxel: brain macromolecular tissue volume in individual subjects. Nature Medicine, 2013 1667-1672.
http://www.nature.com/nm/journal/v19/n12/full/nm.3390.html?WT.ec_id=NM-201312

and the following patent application

>Improved methods for detecting abnormalities in soft tissue using magnetic resonance imaging (MRI),USSN 61/437,587

For more information please contact

>avivmezer@gmail.com



##Contents##


- <a href=#software-requirements>Software Requirements</a>
    - <a href=#required-3rd-party-software>3rd Party software</a>
    - <a href=#matlab-code>Matlab code</a>
- <a href=#mr-scanning->MR Scanning</a>
    - <a href=#spoiled-gradient-echo-scans-spgrflash>Spoiled gradient echo scans (SPGR,FLASH)</a>
    - <a href=#epi-spin-echo-inversion-recovery-scan-b1-mapping>EPI Spin echo inversion recovery scan (B1 mapping)</a>
- <a href=#scanner-dicom-types>Scanner dicom types</a>
- <a href=#data-organization>Data organization</a>
- <a href=#running-mrq>Running mrQ</a>
- <a href=#versions>Versions</a>
- <a href=#parallel-computing>Parallel computing</a>
- <a href=#mrq-analysis-overview>mrQ analysis overview</a>



##Software Requirements##
####Required 3rd party software ####
- MATLAB - http://www.mathworks.com/products/matlab/ 
- ANTS - http://stnava.github.io/ANTs/ 
- FSL - http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/ 
- FreeSrurfer http://surfer.nmr.mgh.harvard.edu/ 
    - *Later versions (>V.1) do not rely on Freesurfer any longer*
- Parallel computing environment (e.g., Sun Grid Engine) 
    - *Not required, but it will make a big difference for running time.*

####Matlab code####
mrQ requires the following openly distributed code repositories:

1. mrQ - https://github.com/mezera/mrQ
2. Vistasoft  - https://github.com/vistalab/vistasoft
3. KNKUtils (from Kendrick Kay) - https://github.com/kendrickkay/knkutils
4. Joëlle Barral's matlab code.  
    - A modified version of this code is distributed within mrQ. 
    - *The original code can be found at: http://www-mrsrl.stanford.edu/~jbarral/t1map.html*

##MR Scanning##
####Spoiled gradient echo scans (SPGR,FLASH)####

1. 2-4 SPGR (not fast SPGR) scans with multiple flip angles recommended (e.g, 4, 10, 20, 30).  
2. All scans should have a single TR (note that a higher TR will increase the SNR).
3. Minimal TE (<=2ms)
4. Save the multi-coil information. To do this on GE scanners, change the scanner default by editing the saveinter cv: saveinter=1.
5. Scan with the same prescan parameters for all SPGR scans. To do this scan the highest SNR image first (flip angle = 10). For the next scan choose manual pre-scan and perform the scan without changing the pre-scan parameters.

####EPI Spin echo inversion recovery scan (B1 mapping)####

Low resolution T1 maps are used to correct for the B1 bias. We will acquire data to fit unbiased T1 maps and correct the bias in the SPGR scans.

The T1 fit is based on Juelle Buarlle’s article: http://onlinelibrary.wiley.com/doi/10.1002/mrm.22497/abstract
http://www-mrsrl.stanford.edu/~jbarral/t1map.html 
A modified version of this code is integrated within the mrQ software. 

1. Scan four SEIR - epi readout scans with four different inversion times (50, 400, 1200, 2400).
2. Each scan needs to be acquired with slab inversion on
GE scanner’s should change the scanner default by editing the a_gzrf0 cv: a_gzrf0=0
3. Use fat suppression. Fat suppression should be spatial-spectral to avoid any slice-selective imperfections. Note: This is the default with GE scanners when slices are less than 4mm thick.

##Scanner dicom types##
The mrQ software was built around GE dicoms. It is possible that different vendors have different conventions in saving dicom information (e.g., header information, data ordering).

We are currently working on making the code compatible with different vendor’s dicom conventions. Please let us know if you experience any issues with reading dicoms when using the software. 

##Data organization##
####Follow these guidelines when organizing your data:####

- Data should be in a single directory - “DATA”.
- Within the DATA directory a dicoms directory is needed.
- Within the dicoms directory each scan should be in a separate directory.
- All SEIR dicoms need to be in the dicom directory.  
- A single dicom is needed for each scan in that scan's directory so that the header can be read by mrQ. 
- All SPGR files need to be in nifti format within the DATA directory.

See http://purl.stanford.edu/qh816pc3429 for an example of directory organization.


##Running mrQ##
####To run mrQ a mrQ structure needs to be created, set and executed.####
For an example of this structure see ‘runScript’ at http://purl.stanford.edu/qh816pc3429

1. Create a structure
    - mrQ=mrQ_Create(path)
2. Set mrQ field 
    - mrQ=mrQ_Set(mrQ,field name,field value)

For a given data set where SEIR scans are organized into 4 folders named  '0005' '0006' '0007' '0008' and SPGR scans are organized into 4 folders named '0009' '0010' '0011' '0012' the following can serve as an example script: 
```matlab
% define the SEIR scans by the session’s 4 characters 
mrQ=mrQ_Set(mrQ,'SEIR',{'0005' '0006' '0007' '0008'})
    
% define the SPGR scans by the session 4 characters
mrQ=mrQ_Set(mrQ,'SPGR',{'0009' '0010' '0011' '0012'})

% make a subject name
mrQ=mrQ_Set(mrQ,'sub','Examp')

% run
mrQ_run(mrQ.name)
```
##Versions##
Version 1 (v.1) is the code to replicate that was used in Nature medicine mezer at. el. 2013 article: https://github.com/mezera/mrQ/tree/v1.0

We recommend you use the most recent, up to date version of mrQ. The most active area of development is in the way the coil sensitivities are calculated. Later versions (>V.1) do not rely on Freesurfer any longer. An article describing those changes is in preparation. 

##Parallel computing##
mrQ takes advantage of parallel computing in three steps within analysis.
1. To calculate transmit inhomogeneity for each voxel.
2. To calculate T1 and M0 to each voxel
3. To calculate the coil gain for different bloc in image space.

mrQ is written to take advantage of the Sun grid parallel computing engine. Each user will need to change the specific calls to the grid according to the parallel computing environment available. One can turn off all those parallel jobs by editing the following setting when creating the mrQ structure:
```matlab
mrQ=mrQ_Set(mrQ,'sungrid’,0);
mrQ=mrQ_Set(mrQ,’proclus’,0);
```
If parallel computing is not available to you please contact us, as we are currently working on a general version of the code that does not rely on parallel computations. 

## T1 fit non linear vs. weighted least square
The most demanding computing is the T1 fit.
to avoid the long .computing (may days on single CPU for all brain  1mm voxel). One can use the  weighted least square as good alternative.

See: Linear least-squares method for unbiased estimation of T1 from SPGR signals. Chang LC, Koay CG, Basser PJ, Pierpaoli C. Magn Reson Med. 2008 Aug;60(2):496-501.

to run it use this setting:
```matlab
mrQ=mrQ_Set(mrQ,'wl’,1);
mrQ=mrQ_Set(mrQ,’lsq’,0);
```
The weighted least square will take few minutes using weighted least square or few hours on a single CPU. 

##mrQ analysis overview##
- mrQ will use the mrQ structure you create and save it to the subject’s directory.
- New directories will be created, including directories for data and quantitative fits.
- Images will be register to each other.
- SEIR-EPI T1 will be computed (low resultion) 
- SPGR T1, M0, B1 maps, and a synthetic T1-weighted image, will be computed.
- T1-weighted and quantitative T1 images will be combined to segment the brain tissue. 
- PD and coil gain will be fit from the M0 image.
- Biophysical model will be applied to calculate VIP and SIR maps.

