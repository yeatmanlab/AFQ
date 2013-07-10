function roi = ANTS_CreateRoiFromMniNifti(roiNifti, warpNifti, outfile)
% Transform a nifti roi saved in MNI space to an individuals brain

% Create the command to call ANTS
cmd = sprintf('WarpImageMultiTransform 3 %s %s -R b0.nii.gz --use-NN %s', roiNifti, warpNifti, outfile)
% excecute it
system(cmd)

% Convert nifti image to .mat roi
roi = dtiRoiFromNifti(outfile,[],prefix(outfile),'mat')