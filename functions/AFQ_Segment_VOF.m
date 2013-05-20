function [vofFG, vofROI1, vofROI2] = AFQ_Segment_VOF(dt, wholebrainFG, varargin)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Define the vertical occipital fasciculus from a wholebrain fiber group.
%
% [vofFG, vofROI1, vofROI2] = AFQ_Segment_VOF(dt, wholebrainFG);
%
% 
% Copyright Jason D Yeatman and Hiromasa Takamura, December 2012


%% Argument checking
if ischar(dt)
    dt = dtiLoadDt6(dt);
end
if ischar(wholebrainFG)
    wholebrainFG = dtiReadFibers(wholebrainFG);
end
% Check if other files were passed in so computations can be skipped
if ~isempty(varargin)
    for ii = 1:length(varargin)
        if isstruct(varargin{ii}) && isfield(varargin{ii},'deformX')...
                && isfield(varargin{ii},'coordLUT')...
                && isfield(varargin{ii},'inMat')
            invDef= varargin{1};
        end
    end
end

% Path to the templates directory
tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');
template = fullfile(tdir,'MNI_EPI.nii.gz');

% Path to the VOF ROIs in MNI space
AFQbase = AFQ_directories;
AFQtemplates = fullfile(AFQbase,'templates');
vof_roi_1 = fullfile(AFQtemplates,'L_VOF_ROI1.nii.gz');
vof_roi_2 = fullfile(AFQtemplates,'L_VOF_ROI2.nii.gz');

% Compute spatial normalization
if ~exist('invDef','var') || isempty(invDef)
    [sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(dt.b0, dt.xformToAcpc, template);
end

% Convert the ROIs to the individual's native space
[~, ~, vofROI1]=dtiCreateRoiFromMniNifti(dt.dataFile, vof_roi_1, invDef);
[~, ~, vofROI2]=dtiCreateRoiFromMniNifti(dt.dataFile, vof_roi_2, invDef);

% Compute the PDD at each point in each ROI
vofROI1_pdd = dtiGetValFromTensors(dt.dt6, vofROI1.coords, inv(dt.xformToAcpc), 'pdd');
vofROI2_pdd = dtiGetValFromTensors(dt.dt6, vofROI2.coords, inv(dt.xformToAcpc), 'pdd');

% Find every coordinate in each VOF ROI where the PDD is in the Z direction
[~, vofROI1_pddZ] = max(abs(vofROI1_pdd),[],2);
vofROI1_pddZ = vofROI1_pddZ == 3;
[~, vofROI2_pddZ] = max(abs(vofROI2_pdd),[],2);
vofROI2_pddZ = vofROI2_pddZ == 3;

% Retain VOF ROI coordinates that have a PDD in the Z direction
vofROI1.coords = vofROI1.coords(vofROI1_pddZ,:);
vofROI2.coords = vofROI2.coords(vofROI2_pddZ,:);

% Intersect the wholebrain fiber group with the ROIs
vofFG = dtiIntersectFibersWithRoi([],'and',[],vofROI1,wholebrainFG);
vofFG = dtiIntersectFibersWithRoi([],'and',[],vofROI2,vofFG);

AFQ_RenderFibers(vofFG,'camera','sagittal')