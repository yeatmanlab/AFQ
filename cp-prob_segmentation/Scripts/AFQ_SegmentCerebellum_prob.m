function afq = AFQ_SegmentCerebellum_prob(afq,afq_fileName,overwriteFiles,sge)
% This function is based on AFQ_SegmentCerebellum:
% 
% September 2021 - Sivan Jossinger:
% Adapted the function for probabilistic tractography.
%
% This function calls cp_AFQ_AddNewFiberGroup and segments the 6 Cerebellar
% Peduncles in probabilistic tractography:
% Superior Cerebellar Peduncle - L&R SCP
% Middle   Cerebellar Peduncle - L&R MCP
% Inferior Cerebellar Peduncle - L&R ICP
%
% inputs:
% afq             - An afq structure output from AFQ_run (see also AFQ_Create)
% overwriteFiles  - 0 or 1. 1 will overwrite to your afq structure.
% sge             - optional parameter. true of false. run on the grid or
% locally. deafault is false- run locally.
%
%  September 2015 - Written by Tal Blecher
%  June 2016      - Editted by Maya Yablonski
%  September 2021 - Editted by Sivan Jossinger 
%
%% Argument checking
if ~exist('overwriteFiles','var') || isempty(overwriteFiles)
    overwriteFiles = [false false];
elseif length(overwriteFiles) == 1;
    if overwriteFiles == 1
        fprintf('\nOVERWRITING ALL PREVIOUS FIBER GROUPS AND ROIS\n');
    end
    overwriteFiles = repmat(overwriteFiles,[1 2]);
elseif length(overwriteFiles) == 2 && overwriteFiles(1) == 0
    fprintf('\nUsing previous CallosumFG but overwriting all segmentations\n');
elseif length(overwriteFiles) > 2
    error('\noverwriteFiles must be a logical vector length 1 or 2\n');
end
if ~exist('sge','var') || isempty(sge)
    % by default run locally not on the grid
    sge = false;
end
if ~exist ('afq_fileName','var') || isempty(afq_fileName)
    d = date;
    afq_fileName = [d '_afq.mat'];
end

% Save cleaning parameters to be able to change back the aqf.params after
% adding the cerebellar tracts
orig_dist = (afq.params.maxDist);
orig_len = (afq.params.maxLen);

% messagebox
mDist = 4; mLen = 1;
message = sprintf('For ICP, cleaning params will be changed to:\nmaxDist: %s \nmaxLen: %s',num2str(mDist), num2str(mLen));
uiwait(msgbox(message, 'Cleaning params changed'));

%% Set up parameters for AFQ_AddNewFiberGroup

% List of ROI names for different cerebellar segments
% ROIs go inferior to superior- to keep order of nodes consistent with
% Travis et al. (2015) HBM paper
% The lateralization is defined with respect to the hemispheres.

roi1Names = {...
    'SCP_L_inferior_prob.nii.gz' 'SCP_R_inferior_prob.nii.gz'...
    'MCP_L_inferior_prob.nii.gz' 'MCP_R_inferior_prob.nii.gz'...
    'ICP_L_inferior_prob.nii.gz' 'ICP_R_inferior_prob.nii.gz' 
    };

roi2Names = {...
    'SCP_R_superior_prob.nii.gz' 'SCP_L_superior_prob.nii.gz'...
    'MCP_R_superior_prob.nii.gz' 'MCP_L_superior_prob.nii.gz'...
    'ICP_L_superior_prob.nii.gz' 'ICP_R_superior_prob.nii.gz' 
    };

% Additional ROI is added to diffrentiate SCP and MCP fibers
roi3Names = {...
    'SCP_L_inter_prob.nii.gz' 'SCP_R_inter_prob.nii.gz'
    };

% Same rois used to define SCP L&R are also used as NOT rois. 
roiNOTNames = {...
    'SCP_L_superior_prob.nii.gz' 'SCP_R_superior_prob.nii.gz'...
    'SCP_L_inter_prob.nii.gz' 'SCP_R_inter_prob.nii.gz'
    };

% Path to the folder containing the template ROIs- or add ROIs to
% mrDiffusion template ROI folder
tdir = uigetdir('','Choose the directory where the Cerebellar ROIs are stored');
for ii = 1:length(roi1Names)
    roi1fullNames{ii} = fullfile(tdir,roi1Names{ii});
    roi2fullNames{ii} = fullfile(tdir,roi2Names{ii});
end
for ii = 1:length(roi3Names)
    roi3fullNames{ii} = fullfile(tdir,roi3Names{ii});
end
for ii = 1:length(roiNOTNames)
    roiNOTfullNames{ii} = fullfile(tdir,roiNOTNames{ii});
end

% List of the cerebellar tracts 
% The lateralization is defined with respect to the cerebellar
% hemispheres.
fgNames = {... 
    'LeftSCP_prob' 'RightSCP_prob'...
    'LeftMCP_prob' 'RightMCP_prob'...
    'LeftICP_prob' 'RightICP_prob'
    };

% Fiber group file to segment
segFgName = 'WholeBrainFG.mat';

%% Segment cerebellum into its different projections and add to afq structure

segNumber=length(afq.fgnames);

SCP_id = find(cellfun('length',regexp(fgNames,'SCP')) == 1);
MCP_id = find(cellfun('length',regexp(fgNames,'MCP')) == 1);
ICP_id = find(cellfun('length',regexp(fgNames,'ICP')) == 1);

for ii = 1:length(fgNames)
    
    if     any(ii == SCP_id)
        afq = cp_AFQ_AddNewFiberGroup4ROI(afq,fgNames{ii},roi1fullNames{ii},roi2fullNames{ii},roi3fullNames{ii},'and',roiNOTfullNames{ii},'not',1,1,0,segFgName,overwriteFiles(2),segNumber+ii);
    
    elseif any(ii == MCP_id)
        afq = cp_AFQ_AddNewFiberGroup3ROI(afq,fgNames{ii},roi1fullNames{ii},roi2fullNames{ii},roiNOTfullNames{ii},'not',1,1,0,segFgName,overwriteFiles(2),segNumber+ii);
    
    elseif any(ii == ICP_id)
        % Change afq cleaning params:
        afq.params.maxDist = mDist;
        afq.params.maxLen = mLen;
        afq = cp_AFQ_AddNewFiberGroup(afq,fgNames{ii},roi1fullNames{ii},roi2fullNames{ii},1,1,0,segFgName,overwriteFiles(2),segNumber+ii);
    else
        error('The fiber group is not a part of the CPs')
    end
    
    fgClipName = [fgNames{ii}(1:end-4) '_clean_D' num2str(afq.params.maxDist) '_L' num2str(afq.params.maxLen) '.mat'];
    ClipFibers(afq,{fgClipName}, roi1Names(ii), roi2Names(ii), tdir, 'and');
end

% Change afq cleaning params to original values:
afq.params.maxDist = orig_dist;
afq.params.maxLen = orig_len;

save(afq_fileName,'afq');

