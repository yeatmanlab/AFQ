function [roi1, roi2] = babyAFQ_LoadROIsReoriented(fgNumber,sub_dir, afq,roiDir)

%%%%%%%%%%%%%%%%%%%%
% This is an infant implementation of AFQ_LoadROIs
% The documentation for the original function in enclosed below
% in contrast to the main babyAFQ_loadROIs, babyAFQ_LoadROIsReoriented
% loads the ROIs in a way that ensure that the axis of the bundles are
% consistent for all fasciles. If this function is combined with
% AFQ_reorientFibers, all fascicles will be oriented posterior -> anterior,
% inferior -> superior and left -> right

% In brief, the implemented infant-specfic changes are:
% 1) additional ROIs were added, so that ROIs for pAF and mdLF can also be
% loaded with this function

% InfantAFQ was written by Mareike Grotheer in 2020

%%%%%%%%%%%%%%%%%%%%%

% Load the two ROIs used to define a fiber group
%
% [roi1 roi2] = AFQ_LoadROIs(fgNumber,sub_dir. [afq])
%
% Inputs:
% fgNumber = The number of the fiber group. For example 1 is left ATR and
%            20 is right arcuate.
% sub_dir  = A path to the subjects directory
% afq      - An afq structure. This is optional but is needed if it is a
%            fiber group not in the mori groups
%
% Outputs:
% roi1     = First roi defining the tract
% roi2     = Second roi defining the tract
%
% Written by Jason D. Yeatman 1/31/2012

% If no afq structure was passed in it is fine for known ROIs
if ~exist('afq','var') || isempty(afq) && fgNumber <25
    roi1Names={'ATR_roi2_L', 'ATR_roi2_R','CST_roi1_L','CST_roi1_R','CGC_roi1_L','CGC_roi1_R'...
        'HCC_roi1_L','HCC_roi1_R','FP_L','FA_L','IFO_roi1_L','IFO_roi1_R','ILF_roi1_L','ILF_roi1_R'...
        'SLF_roi2_L','SLF_roi2_R','UNC_roi1_L','UNC_roi1_R','SLFt_roi2_L','SLFt_roi2_R'...
        'MdLF_roi1_L','MdLF_roi1_R','SLFt_roi2_L','SLFt_roi2_R'};
    roi2Names={'ATR_roi1_L', 'ATR_roi1_R','CST_roi2_L','CST_roi2_R','CGC_roi2_L','CGC_roi2_R'...
        'HCC_roi2_L','HCC_roi2_R','FP_R','FA_R','IFO_roi2_L','IFO_roi2_R','ILF_roi2_L','ILF_roi2_R'...
        'SLF_roi1_L','SLF_roi1_R','UNC_roi2_L','UNC_roi2_R','SLF_roi1_L','SLF_roi1_R'...
        'ILF_roi2_L','ILF_roi2_R','LH_Parietal', 'RH_Parietal'};
    
    roi1path = fullfile(roiDir,[roi1Names{fgNumber} '.mat']);
    roi2path = fullfile(roiDir,[roi2Names{fgNumber} '.mat']);
    if exist(roi1path,'file') && exist(roi2path,'file')
        roi1=dtiReadRoi(roi1path);
        roi2=dtiReadRoi(roi2path);
    else
        error('ROI does not exist for %s',sub_dir)
    end
else
    % It doesn't matter which subject number because all we need is the roi
    % name
    try
    [~,roi1Name] = fileparts(AFQ_get(afq,'roi1',fgNumber,1));
    [~,roi2Name] = fileparts(AFQ_get(afq,'roi2',fgNumber,1));

    roi1path = fullfile(roiDir,[roi1Name '.mat']);
    roi2path = fullfile(roiDir,[roi2Name '.mat']);
    if exist(roi1path,'file') && exist(roi2path,'file')
        roi1=dtiReadRoi(roi1path);
        roi2=dtiReadRoi(roi2path);
    else
        error('ROI does not exist for %s',sub_dir)
    end
    catch
        roi1=[];
        roi2=[];
    end
end

return