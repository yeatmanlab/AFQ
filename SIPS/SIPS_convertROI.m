function SIPS_convertROI(fsFolder, anatomyFile)

% This code aims to generate ROI files requires for human SIPS identification.
% Procedure: 
% (1) Convert freesurfer-based .mat ROI into nifti format.
% (2) Merge Precuneus and Superior Parietal ROI for SIPS identification
%
% Requirement:
% (1) Free-surfer segmentation: 
% e.g) recon-all -i t1.nii.gz -subjid HCP_001 -all 
% (2) Convert Freesurfer ROI into .mat format (use fs_roiFromAllLabels in vistasoft repository)
%
% Dependency:
% vistasoft: https://github.com/vistalab/vistasoft
% INPUT:
% fsFolder: The full path to folder with freesurfer ROI in .mat format
% anatomyFile: The full path to anatomy file (T1-weighted) for a reference. Anatomy file must be coregistred with diffusion data. 
% 
% If you use this code for your own study, please cite following article as a reference:
% Uesaki, M., Takemura, H. & Ashida, H. (2017) Computational neuroanatomy
% of human stratum proprium of interparietal sulcus. Brain Structure and
% Function, in press.
%
% Documentation: https://github.com/vistalab/vistasoft/wiki/Identify-human-Stratum-Proprium-of-Interparietal-Sulcus
% 
% (C) Hiromasa Takemura, CiNet NICT HHS 2017

parietalROIs = {'1025_ctx-lh-precuneus.mat','1029_ctx-lh-superiorparietal.mat',...
    '1031_ctx-lh-supramarginal.mat','2025_ctx-rh-precuneus.mat','2029_ctx-rh-superiorparietal.mat',...
    '2031_ctx-rh-supramarginal.mat'};

cd(fsFolder);

% Convert ROI into nifti
for i = 1:length(parietalROIs)
dtiRoiNiftiFromMat(parietalROIs{i},anatomyFile);
end

% Load rois
lh_precunueus = niftiRead('1025_ctx-lh-precuneus.nii.gz');
lh_spl = niftiRead('1029_ctx-lh-superiorparietal.nii.gz');
lh_smg = niftiRead('1031_ctx-lh-supramarginal.nii.gz');
rh_precunueus = niftiRead('2025_ctx-rh-precuneus.nii.gz');
rh_spl = niftiRead('2029_ctx-rh-superiorparietal.nii.gz');
rh_smg = niftiRead('2031_ctx-rh-supramarginal.nii.gz');

% Merge SPL and Precuneus, for left hemisphere
lh_new2 = lh_smg; 
lh_new2.data = zeros(size(lh_smg.data));
lh_new2.data(lh_precunueus.data == 1) = 1;
lh_new2.data(lh_spl.data == 1) = 1;
lh_new2.fname = 'lh_precuneus_spl.nii.gz';
niftiWrite(lh_new2);

% Merge SPL and Precuneus, for right hemisphere
rh_new2 = rh_smg; 
rh_new2.data = zeros(size(rh_smg.data));
rh_new2.data(rh_precunueus.data == 1) = 1;
rh_new2.data(rh_spl.data == 1) = 1;
rh_new2.fname = 'rh_precuneus_spl.nii.gz';
niftiWrite(rh_new2);
