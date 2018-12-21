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




