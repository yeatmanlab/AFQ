function [L_FG R_FG L_roi_2 L_roi_3 R_roi_2 R_roi_3] = AFQ_Segment_PostArcuate(dt, wholebrainFG, sub_dir, varargin)
% Define the posterior segment of the arcuate from a wholebrain fiber group
%
% [L_FG R_FG L_roi_2 L_roi_3 R_roi_2 R_roi_3] = AFQ_Segment_PostArcuate(dt, wholebrainFG, sub_dir, varargin)
%

%% Argument checking
if ~exist('dt','var') || isempty(dt)
    error('dt6 file is needed');
elseif ischar(dt)
    dt = dtiLoadDt6(dt);
end
if ~exist('wholebrainFG','var') || isempty(wholebrainFG)
    fprintf('Wholebrain fiber group was not supplied. Tracking fibers');
    wholebrainFG = AFQ_WholebrainTractography(dt);
elseif ischar(wholebrainFG) && exist(wholebrainFG,'file')
    wholebrainFG = dtiReadFibers(wholebrainFG);
elseif ischar(wholebrainFG) && ~exist(wholebrainFG,'file')
    fprintf('Wholebrain fiber group was not found. Tracking fibers');
    wholebrainFG = AFQ_WholebrainTractography(dt);
end
if ~exist('sub_dir','var') || isempty(sub_dir)
    sub_dir = fileparts(dt.dataFile);
end
% Check if inverse deformation was passed in so computations can be skipped
if ~isempty(varargin)
    for ii = 1:length(varargin)
        if isstruct(varargin{ii}) && isfield(varargin{ii},'deformX')...
                && isfield(varargin{ii},'coordLUT')...
                && isfield(varargin{ii},'inMat')
            invDef= varargin{1};
        end
    end
end

%% Create ROIs and segment the posterior segment of the arcuate
% Path to the templates directory
tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');
template = fullfile(tdir,'MNI_EPI.nii.gz');

% Path to the VOF ROIs in MNI space
AFQbase = AFQ_directories;
L_roi_2 = fullfile(AFQbase,'templates', 'L_Arcuate_PostSegment.mat');
R_roi_2 = fullfile(AFQbase,'templates', 'R_Arcuate_PostSegment.mat');
L_roi_3 = fullfile(AFQbase,'templates', 'L_LateralParietal.mat');
R_roi_3 = fullfile(AFQbase,'templates', 'R_LateralParietal.mat');

% Compute spatial normalization
if ~exist('invDef','var') || isempty(invDef)
    [sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(dt.b0, dt.xformToAcpc, template);
end

% Load up template ROIs in MNI space
L_roi_2 = dtiReadRoi(L_roi_2);
R_roi_2 = dtiReadRoi(R_roi_2);
L_roi_3 = dtiReadRoi(L_roi_3);
R_roi_3 = dtiReadRoi(R_roi_3);
% Transform ROI coords to the subjects native space
L_roi_2.coords = mrAnatXformCoords(invDef,L_roi_2.coords);
R_roi_2.coords = mrAnatXformCoords(invDef,R_roi_2.coords);
L_roi_3.coords = mrAnatXformCoords(invDef,L_roi_3.coords);
R_roi_3.coords = mrAnatXformCoords(invDef,R_roi_3.coords);
% Save the rois
dtiWriteRoi(L_roi_2,fullfile(sub_dir,'ROIs',L_roi_2.name));
dtiWriteRoi(R_roi_2,fullfile(sub_dir,'ROIs',R_roi_2.name));
dtiWriteRoi(L_roi_3,fullfile(sub_dir,'ROIs',L_roi_3.name));
dtiWriteRoi(R_roi_3,fullfile(sub_dir,'ROIs',R_roi_3.name));
% Load up predefined ROIs for the Arcuate fasciculus. These are computed in
% AFQ_SegmentFiberGroups
[L_roi_not, L_roi_1] = AFQ_LoadROIs(19,sub_dir);
[R_roi_not, R_roi_1] = AFQ_LoadROIs(20,sub_dir);

% Intersect the wholebrain fiber group with the ROIs to retain only fibers
% that correspond to the posterior segment of the arcuate
L_FG = dtiIntersectFibersWithRoi([],'and',[],L_roi_1,wholebrainFG);
L_FG = dtiIntersectFibersWithRoi([],'and',[],L_roi_2,L_FG);
L_FG = dtiIntersectFibersWithRoi([],'and',[],L_roi_3,L_FG);
L_FG = dtiIntersectFibersWithRoi([],'not',[],L_roi_not,L_FG);
R_FG = dtiIntersectFibersWithRoi([],'and',[],R_roi_1,wholebrainFG);
R_FG = dtiIntersectFibersWithRoi([],'and',[],R_roi_2,R_FG);
R_FG = dtiIntersectFibersWithRoi([],'and',[],R_roi_3,R_FG);
R_FG = dtiIntersectFibersWithRoi([],'not',[],R_roi_not,R_FG);

% Name the fiber groups
L_FG.name = 'L_Arcuate_Posterior';
R_FG.name = 'R_Arcuate_Posterior';

% Save the segmented fiber groups
dtiWriteFiberGroup(L_FG,fullfile(sub_dir,'fibers',L_FG.name));
dtiWriteFiberGroup(R_FG,fullfile(sub_dir,'fibers',R_FG.name))

% Show the fiber group if desired
if sum(strcmpi('showfibers', varargin)) > 0
    AFQ_RenderFibers(L_FG,'color',[1 0 0],'numfibers',300);
    b0 = readFileNifti(dt.files.b0);
    AFQ_AddImageTo3dPlot(b0,[-30 0 0]);
end
return