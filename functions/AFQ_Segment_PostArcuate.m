function [L_FG, R_FG, L_roi_1, L_roi_2, R_roi_1, R_roi_2] = AFQ_Segment_PostArcuate(dt, wholebrainFG, varargin)
% Define the posterior segment of the arcuate from a wholebrain fiber group
%
% [L_FG, R_FG, L_roi_1, L_roi_2, R_roi_1, R_roi_2] = AFQ_Segment_PostArcuate(dt, wholebrainFG, varargin)
%
% This function will define the left and right posterior segment of the
% arcuate fasciculus. 
%
% Arguments:
% 'showfibers' - Locigal. Whether or not to render figures of the fibers
% 'saveFiles'  - Logical. Whether or not to save the fibers and ROIs
%
% Copyright Jason D. Yeatman December 2012
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

if ~isempty(varargin)
    % Check if inverse deformation was passed in so computations can be skipped
    for ii = 1:length(varargin)
        if isstruct(varargin{ii}) && isfield(varargin{ii},'deformX')...
                && isfield(varargin{ii},'coordLUT')...
                && isfield(varargin{ii},'inMat')
            invDef= varargin{ii};
        end
    end
    % Check if rois and fibers should be saved
    for ii = 1:length(varargin)
        if strcmpi('saveFiles',varargin{ii}) && length(varargin) > ii
            saveFiles = varargin{ii+1};
        end
    end
end

% Save figures by default
if ~exist('saveFiles','var') || isempty(saveFiles)
    saveFiles = 1;
end

%% Create ROIs and segment the posterior segment of the arcuate
% Path to the templates directory
tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');
template = fullfile(tdir,'MNI_EPI.nii.gz');

% Path to the VOF ROIs in MNI space
AFQbase = AFQ_directories;
% L_roi_2 = fullfile(AFQbase,'templates', 'L_Arcuate_PostSegment.nii.gz');
% R_roi_2 = fullfile(AFQbase,'templates', 'R_Arcuate_PostSegment.nii.gz');
L_roi_2 = fullfile(AFQbase,'templates', 'L_Parietal.nii.gz');
R_roi_2 = fullfile(AFQbase,'templates', 'R_Parietal.nii.gz');
% L_roi_2 = fullfile(AFQbase,'templates', 'L_LateralParietal.nii.gz');
% R_roi_2 = fullfile(AFQbase,'templates', 'R_LateralParietal.nii.gz');

% Compute spatial normalization
if ~exist('invDef','var') || isempty(invDef)
    [sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(dt.b0, dt.xformToAcpc, template);
end

% Load up template ROIs in MNI space and transform them to the subjects
% native space.  dtiCreateRoiFromMniNifti also saves the ROIs
[~, ~, L_roi_2]=dtiCreateRoiFromMniNifti(dt.dataFile, L_roi_2, invDef, 1);
[~, ~, R_roi_2]=dtiCreateRoiFromMniNifti(dt.dataFile, R_roi_2, invDef, 1);
% Name the ROIs
L_roi_2.name = 'L_Parietal';
R_roi_2.name = 'R_Parietal';

% Load up predefined ROIs for the Arcuate fasciculus. These are computed in
% AFQ_SegmentFiberGroups
[L_roi_not, L_roi_1] = AFQ_LoadROIs(19,sub_dir);
[R_roi_not, R_roi_1] = AFQ_LoadROIs(20,sub_dir);

% To eliminate fibers that head too far anterior we will create a full
% plane at the y coordinate of the anterior slf roi (around the central
% sulcus) and remove fibers that intersect this ROI
L_roi_not = dtiRoiMakePlane([-60,mean(L_roi_not.coords(:,2)),-60; 60,mean(L_roi_not.coords(:,2)),80]);
R_roi_not = dtiRoiMakePlane([-60,mean(R_roi_not.coords(:,2)),-60; 60,mean(R_roi_not.coords(:,2)),80]);

% Intersect the wholebrain fiber group with the ROIs to retain only fibers
% that correspond to the posterior segment of the arcuate
L_FG = dtiIntersectFibersWithRoi([],'and',[],L_roi_1,wholebrainFG);

% We want fibers that are traveling superior to inferior as they pass
% through the first ROI.  To do this we will restrict the ROI to voxels
% with superior-inferior orientation and repeate this process for planes
% directly above and below the ROI
for ii = -3:3
    roi_tmp = L_roi_1;
    % Shift the roi
    roi_tmp.coords(:,3) = roi_tmp.coords(:,3) + ii;
    % Compute the PDD at each point in each ROI
    roi_pdd = dtiGetValFromTensors(dt.dt6, roi_tmp.coords, inv(dt.xformToAcpc), 'pdd');
    % Find every coordinate in each VOF ROI where the PDD is in the Z direction
    [~, roi_pddZ] = max(abs(roi_pdd),[],2);
    roi_pddZ = roi_pddZ == 3;
    % Retain VOF ROI coordinates that have a PDD in the Z direction
    roi_tmp.coords = roi_tmp.coords(roi_pddZ,:);
    L_FG = dtiIntersectFibersWithRoi([],'and',[],roi_tmp,L_FG);
end
%L_FG = dtiIntersectFibersWithRoi([],'and',[],L_roi_2,L_FG);
% Make sure endpoints are in the parietal lobe
L_FG = dtiIntersectFibersWithRoi([],'endpoints',4,L_roi_2,L_FG);
% And that fibers do not head to far anterior
L_FG = dtiIntersectFibersWithRoi([],'not',[],L_roi_not,L_FG);

% Same for the right hemisphere fiber group
R_FG = dtiIntersectFibersWithRoi([],'and',[],R_roi_1,wholebrainFG);
for ii = -3:3
    roi_tmp = R_roi_1;
    % Shift the roi
    roi_tmp.coords(:,3) = roi_tmp.coords(:,3) + ii;
    % Compute the PDD at each point in each ROI
    roi_pdd = dtiGetValFromTensors(dt.dt6, roi_tmp.coords, inv(dt.xformToAcpc), 'pdd');
    % Find every coordinate in each VOF ROI where the PDD is in the Z direction
    [~, roi_pddZ] = max(abs(roi_pdd),[],2);
    roi_pddZ = roi_pddZ == 3;
    % Retain VOF ROI coordinates that have a PDD in the Z direction
    roi_tmp.coords = roi_tmp.coords(roi_pddZ,:);
    R_FG = dtiIntersectFibersWithRoi([],'and',[],roi_tmp,R_FG);
end
%R_FG = dtiIntersectFibersWithRoi([],'and',[],R_roi_2,R_FG);
% Make sure endpoints are in the parietal lobe
R_FG = dtiIntersectFibersWithRoi([],'endpoints',4,R_roi_2,R_FG);
% And that fibers do not head to far anterior
R_FG = dtiIntersectFibersWithRoi([],'not',[],R_roi_not,R_FG);

% Name the fiber groups
L_FG.name = 'L_Arcuate_Posterior';
R_FG.name = 'R_Arcuate_Posterior';

% Save the segmented fiber groups and ROIs
if saveFiles == 1
    dtiWriteFiberGroup(L_FG,fullfile(sub_dir,'fibers',L_FG.name));
    dtiWriteFiberGroup(R_FG,fullfile(sub_dir,'fibers',R_FG.name));
    dtiWriteRoi(L_roi_2,fullfile(sub_dir,'ROIs',L_roi_2.name));
    dtiWriteRoi(R_roi_2,fullfile(sub_dir,'ROIs',R_roi_2.name));
end

% Show the fiber group if desired
if sum(strcmpi('showfibers', varargin)) > 0
    AFQ_RenderFibers(L_FG,'color',[0 0 1],'numfibers',300);
    b0 = readFileNifti(dt.files.b0);
    AFQ_AddImageTo3dPlot(b0,[-30 0 0]);
    AFQ_RenderRoi(L_roi_1);
end

return