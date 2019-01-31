function afq = AFQ_SegmentCerebellum(afq,afq_fileName,overwriteFiles,sge)

%% This function calls cp_AFQ_AddNewFiberGroup/cp_AFQ_AddNewFiberGroup3ROI and segments the cerebellar peduncles:
%  ICP - left/right inferior cerebellar peduncle
%  SCP - left/right superior cerebellar peduncle
%  MCP - middle cerebellar peduncle

%  Written by Tal Blecher - September 2015
%  Adapted by Maya Yablonski - June 2016 
%  Edited for publication - January 2019

%% Input variables:
%  afq             - An afq structure output from AFQ_run (see also AFQ_Create)
%  overwriteFiles  - Options: 0 = do not overwrite existing afq structure (default)
%                             1 = overwrite existing afq structure
%  sge             - Options: false = run locally (default)
%                             true = run on the grid

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
    sge = false;
end
if ~exist ('afq_fileName','var') || isempty(afq_fileName)
    d = date;
    afq_fileName = [d '_afq.mat'];
end

%% Cleaning parameters for segmented fiber groups

%  Save default afq cleaning parameters so that they can be restored after cerebellar segmentation
orig_dist = (afq.params.maxDist);
orig_len = (afq.params.maxLen);

%  Change cleaning parameter for cerebellar segmentation and display messagebox
mDist = 4; 
mLen = 1;

message = sprintf('For Cerebellar fibers, cleaning params will be changed to:\nmaxDist: %s \nmaxLen: %s',num2str(mDist), num2str(mLen));
uiwait(msgbox(message, 'Cleaning params changed'));

afq.params.maxDist = mDist;
afq.params.maxLen = mLen;

%% Set up parameters for running modified AFQ_AddNewFiberGroup functions

% List of ROI names for different cerebellar segments
% Note: ROIs go from bottom to top and from left to right to keep
% the order of nodes consistent with Travis et al. (2015, HBM)
roi1Names = {'ICP_L_low.nii.gz' 'ICP_R_low.nii.gz' ...
             'SCP_L_low.nii.gz' 'SCP_R_low.nii.gz' 'MCP_L.nii.gz'};
roi2Names = {'ICP_L_high.nii.gz' 'ICP_R_high.nii.gz' ...
             'SCP_L_high.nii.gz' 'SCP_R_high.nii.gz' 'MCP_R.nii.gz'};

% Additional ROI to remove MCP streamlines crossing the midline in the pons
MCP_NOT_roiFileName = 'midline_wedgebig.nii.gz';

% Specify path to the cerebellar ROIs
tdir = uigetdir('','Choose the directory where the Cerebellar ROIs are stored');
for ii = 1:length(roi1Names)
    roi1fullNames{ii} = fullfile(tdir,roi1Names{ii});
    roi2fullNames{ii} = fullfile(tdir,roi2Names{ii});
end

MCP_NOT_roiFileName = fullfile(tdir, MCP_NOT_roiFileName);

% Alternatively, one can add the cerebellar ROIs to the mrDiffusion template ROI folder 
% and uncomment the following line:
% tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');

% List of the cerebellar peduncle names
fgNames = {'LeftICP.mat' 'RightICP.mat' ...
           'LeftSCP.mat' 'RightSCP.mat' 'MCP.mat'};

% Specify whole-brain tractogram used to segment the cerebellar peduncles
segFgName = 'WholeBrainFG.mat';

%% Segment cerebellar peduncles and add them to afq structure

% Call modified AFQ_AddNewFiberGroup functions and clip fibers
segNumber=length(afq.fgnames);

for ii = 1:length(fgNames)
    MCP_id = find(cellfun('length',regexp(fgNames,'MCP')) == 1);
    if ii == MCP_id
        afq = cp_AFQ_AddNewFiberGroup3ROI(afq,fgNames{ii},roi1fullNames{ii},roi2fullNames{ii},MCP_NOT_roiFileName,'not',1,1,0,segFgName,overwriteFiles(2),segNumber+ii);
    else
        afq = cp_AFQ_AddNewFiberGroup(afq,fgNames{ii},roi1fullNames{ii},roi2fullNames{ii},1,1,0,segFgName,overwriteFiles(2),segNumber+ii);
    end
    fgClipName = [fgNames{ii}(1:end-4) '_clean_D' num2str(mDist) '_L' num2str(mLen) '.mat'];
    ClipFibers(afq,{fgClipName}, roi1Names(ii), roi2Names(ii), tdir, 'and');
end

afq_newName = [afq_fileName(1:end-4) '_cerebellar.mat'];

% Change afq cleaning parameters back to default values:
afq.params.maxDist = orig_dist;
afq.params.maxLen = orig_len;

save(afq_newName,'afq');