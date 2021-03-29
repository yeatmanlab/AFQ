function img=  babyAFQ_PrepareVofROI(fatDir,sessid,runName,ROIdir,asegFile, fivettFile)
%This function prepares the volumetric VOF ROI by intersecting it with the
%Gray-matter/white-matter interface generated from mrtrix
babyAFQ_DtiRoi2Nii(fatDir, sessid, runName, []);
runDir=fullfile(fatDir,sessid,runName);
cmd_str=['mri_convert ' fullfile(ROIdir,'VOF_box_L.nii.gz') ' ',...
    fullfile(ROIdir,'VOF_box_L_resliced.nii.gz') ' --reslice_like ',...
    fullfile(asegFile)];
system(cmd_str);

cmd_str=['mri_convert ' fullfile(ROIdir,'VOF_box_R.nii.gz') ' ',...
    fullfile(ROIdir,'VOF_box_R_resliced.nii.gz') ' --reslice_like ',...
    fullfile(asegFile)];
system(cmd_str);

cmdstr=['5tt2gmwmi -mask_in ' fullfile(ROIdir,'VOF_box_L_resliced.nii.gz') ' ',...
    fullfile(fivettFile) ' ',...
    fullfile(ROIdir,'VOF_gmwmi_L.nii.gz') ' -force'];
system(cmdstr);

cmdstr=['5tt2gmwmi -mask_in ' fullfile(ROIdir,'VOF_box_R_resliced.nii.gz') ' ',...
    fullfile(fivettFile) ' ',...
    fullfile(ROIdir,'VOF_gmwmi_R.nii.gz') ' -force'];
system(cmdstr);

img=niftiRead(fullfile(ROIdir,'VOF_gmwmi_L.nii.gz'));
load(fullfile(ROIdir,'VOF_box_L.mat'));

imgCoords = find(ceil(img.data));
[I,J,K] = ind2sub(img.dim, imgCoords);
roi.coords  = mrAnatXformCoords(img.qto_xyz, [I,J,K]);

outname='VOF_gmwmi_L.mat';
save(fullfile(ROIdir,outname),'roi','versionNum','coordinateSpace');

img=niftiRead(fullfile(ROIdir,'VOF_gmwmi_R.nii.gz'));
load(fullfile(ROIdir,'VOF_box_R.mat'));

imgCoords = find(ceil(img.data));
[I,J,K] = ind2sub(img.dim, imgCoords);
roi.coords  = mrAnatXformCoords(img.qto_xyz, [I,J,K]);

outname='VOF_gmwmi_R.mat';
save(fullfile(ROIdir,outname),'roi','versionNum','coordinateSpace');
end