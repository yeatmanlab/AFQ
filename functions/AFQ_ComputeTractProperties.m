function [fa md rd ad cl TractProfile] = AFQ_ComputeTractProperties(fg_classified,dt,numberOfNodes,clip2rois,subDir)
% Compute diffusion properties along the trajectory of the fiber groups
%
% [fa md rd ad cl TractProfile] = AFQ_ComputeTractProperties(fg_classified,dt,[numberOfNodes=30],[clip2rois=1],[subDir])
%
% Input arguments:
%  fg_classified  = Fiber group classified into 20 subgroups that
%                   correspond to each of the 20 tracts output by
%                   AFQ_SegmentFiberGroups
%
%  dt             = dt6 structure or an image file.
%
%  numberOfNodes  = The number of nodes to sample each fiber to.
%
%  clip2rois      = if clip2rois is set to 1 (default) then only the
%                   central portion of the fiber group spanning between
%                   the 2 defining ROIs is analyzed.
%                   Otherwise the full trajectory is analyzed.
%
%  subDir         = The directory containing the subjects dt6 file. This is
%                   needed to find the subject's ROIs if clip2rois==1. By
%                   defualt the subject's directory will be searched for
%                   based on the path given in the dt6 file (dt.dataFile).
%
% Outputs:
% fa               = vector of Fractional anisotropy along fiber core.
% md               = vector of Mean Diffusivity values along fiber core.
% rd               = vector of Radial Diffusivity values along fiber core.
% ad               = vector of Axial Diffusivity values along fiber core.
% cl               = vector of linearity values along fiber core.
% SuperFibersGroup = Fiber group core and covariance at each node
%
%  Example:
%
%   AFQdata = '/home/jyeatman/matlab/svn/vistadata/AFQ'
%   dt = dtiLoadDt6(fullfile(AFQdata,'subj2','dt6.mat'))
%   fg_classified = AFQ_SegmentFiberGroups(dt)
%   numberOfNodes = 30; clip2rois = 1; subDir = fullfile(AFQdata,'subj2');
%   [fa md rd ad cl SuperFibersGroup] =
%   AFQ_ComputeTractProperties(fg_classified,dt,numberOfNodes,clip2rois,subDir)
%
% Copyright Stanford Vista Team 2011
% Written by Jason D. Yeatman

% check the format of the fiber group.
if isstruct(fg_classified) && length(fg_classified) > 1
    fg_classified = dtiFgArrayToFiberGroup(fg_classified);
end
% if number of nodes doesn't exist then set to 30
if ~exist('numberOfNodes','var') || isempty(numberOfNodes)
    numberOfNodes=30;
end
if ~exist('clip2rois','var') || isempty(clip2rois)
    clip2rois=1;
end
% check if the dt input is a dt6 structure or a nifti image
if isfield(dt,'qto_xyz')
    valname = 'image';
else
    valname = 'dt';
end
% if no subject directory was defined then define based on the dt6 file
if ~exist('subDir','var') || isempty(subDir) && strcmp(valname, 'dt')
    subDir = fileparts(dt.dataFile);
end
% Create Tract Profile structure
TractProfile = AFQ_CreateTractProfile;
% loop over the fiber groups
for jj=1:length(fg_classified.subgroupNames)
    % There are 20 fiber groups saved with the same structure. We will loop
    % over these fiber groups, allocate them to a tmp variable calculate
    % what we need and move on
    fgtmp=dtiNewFiberGroup(fg_classified.subgroupNames(jj).subgroupName);
    fgtmp.fibers=fg_classified.fibers(fg_classified.subgroup==jj);
    % clip the fiber group to the portion spanning between the two ROIs if
    % desired
    if clip2rois==1
        [roi1 roi2] = AFQ_LoadROIs(jj,subDir);
    else
        roi1=[]; roi2=[];
    end
    % Some subjects may be missing certain fiber groups due to lesions,
    % or noise in the data.  If a fiber group is empty we will ad NaNs to
    % the database
    if isempty(fgtmp.fibers)
        fa(:, jj)=nan;md(:, jj)=nan;rd(:, jj)=nan;ad(:, jj)=nan;cl(:, jj)=nan;
        TractProfile(jj) = AFQ_CreateTractProfile('name',fgtmp.name);
        continue
    end
    % Here is where we will create define the fiber group core and calculate
    % stats based on the method described in Yeatman et al 2011. Sometimes
    % there is a problem such as only having 1 fiber in the fiber group.
    % To avoid breaking the whole loop we use the try catch loop
    try
        [fa(:, jj),md(:, jj),rd(:, jj),ad(:, jj),cl(:, jj),SuperFibersGroup(jj)]=dtiComputeDiffusionPropertiesAlongFG(fgtmp, dt, roi1, roi2, numberOfNodes);
        % Put mean fiber into Tract Profile structure
        TractProfile(jj) = AFQ_CreateTractProfile('name',fgtmp.name,'superfiber',SuperFibersGroup(jj));
    catch
        fa(:, jj)=nan;md(:, jj)=nan;rd(:, jj)=nan;ad(:, jj)=nan;cl(:, jj)=nan;
        TractProfile(jj) = AFQ_CreateTractProfile('name',fgtmp.name);
    end
end


return