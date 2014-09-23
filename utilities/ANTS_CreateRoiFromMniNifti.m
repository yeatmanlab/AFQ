function roi = ANTS_CreateRoiFromMniNifti(roiNifti, warpNifti,refNifti, outfile)
% Transform a nifti roi saved in MNI space to an individuals brain

% Find the affine associated with the warp
f = strfind(warpNifti,'Inverse')-1;
affineXform = [warpNifti(1:f) 'Affine.txt'];

% Create the command to call ANTS
if exist('refNifti','var') && ~isempty(refNifti)
    % If a reference image was provided then reslice the roi to match this
    cmd = sprintf('WarpImageMultiTransform 3 %s %s -R %s --use-NN -i %s %s', roiNifti, outfile, refNifti, affineXform, warpNifti);
else
    cmd = sprintf('WarpImageMultiTransform 3 %s %s --use-NN -i %s %s', roiNifti, outfile, affineXform, warpNifti);
end
% excecute it
system(cmd);

% Convert nifti image to .mat roi
roi = dtiRoiFromNifti(outfile,[],prefix(prefix(outfile)),'mat');

% Clean the ROI
roi = dtiRoiClean(roi,3,{'fillHoles'});

% Save the ROI
dtiWriteRoi(roi,prefix(prefix(outfile)));