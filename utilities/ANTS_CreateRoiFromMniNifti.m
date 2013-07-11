function roi = ANTS_CreateRoiFromMniNifti(roiNifti, warpNifti,refNifti, outfile)
% Transform a nifti roi saved in MNI space to an individuals brain

% Create the command to call ANTS
if exist('refNifti','var') && ~isempty(refNifti)
    % If a reference image was provided then reslice the roi to match this
    cmd = sprintf('WarpImageMultiTransform 3 %s %s -R %s --use-NN %s', roiNifti, outfile, refNifti, warpNifti);
else
    cmd = sprintf('WarpImageMultiTransform 3 %s %s --use-NN %s', roiNifti, outfile, warpNifti);
end
% excecute it
system(cmd);

% Convert nifti image to .mat roi
roi = dtiRoiFromNifti(outfile,[],prefix(prefix(outfile)),'mat');

% Clean the ROI
roi = dtiRoiClean(roi,3,{'dilate' 'fillHoles'});

% Save the ROI
dtiWriteRoi(roi,prefix(prefix(outfile)));