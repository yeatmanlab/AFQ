function [volRoi] = AFQ_meshDrawRoiVolFill(msh, voldata, fill_range, fh)
%
%
%
% Example:
%
% voldata = 'functionalOverlay-20-Dec-2018.nii.gz';
% t1class = 't1_class.nii.gz';
% im = niftiRead(t1class);
% msh = AFQ_meshCreate(t1class);
% fill_range = [2 inf];
% volRoi = AFQ_meshDrawRoiVolFill(msh, voldata, fill_range);

if ~exist('fh', 'var')
    fh = [];
end

voldata = niftiRead(voldata); voldata.fname = 'voldata.nii.gz';
msh = AFQ_meshColor(msh,'overlay',voldata, 'thresh',fill_range);
[coords, indices, bin] = AFQ_meshDrawRoi(msh, 0, fh);
imcoords = ceil(mrAnatXformCoords(voldata.qto_ijk, coords));
v = voldata.data > fill_range(1) & voldata.data<fill_range(2);
volRoi = imfill(v,imcoords');


