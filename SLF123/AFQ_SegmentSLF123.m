function afq = AFQ_SegmentSLF123(afq,afq_fileName,overwriteFiles,sge)

% August 2023 - Romi Sagi:
% This function calls SLF123_AFQ_AddNewFiberGroup3ROI and segments the 3
% Superior Longitudinal Fasciculus branches - SLF-I,-II, and -III.

% Inputs:
% afq             - An afq structure output from AFQ_run (see also AFQ_Create)
% afq_fileName    - A string containing the name under which the new AFQ structure will be saved
% overwriteFiles  - Logical indicating whether ot not to overwrite to your afq structure
% sge             - Logical defining whether to run on the grid or locally. Deafault is false- 
%                   run locally (optional parameter)

%% Argument checking
if ~exist('overwriteFiles','var') || isempty(overwriteFiles)
    overwriteFiles = false;
elseif length(overwriteFiles) == 1
    if overwriteFiles == 1
        fprintf('\nOVERWRITING ALL PREVIOUS FIBER GROUPS AND ROIS\n');
    end
elseif length(overwriteFiles) > 1
    error('\noverwriteFiles must be a logical vector length 1\n');
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
% adding the SLF branches
orig_dist = (afq.params.maxDist);
orig_len = (afq.params.maxLen);

% messagebox
afq.params.maxDist = 4; afq.params.maxLen = 1;
message = sprintf('For SLF123, cleaning params will be changed to:\nmaxDist: %s \nmaxLen: %s',num2str(mDist), num2str(mLen));
uiwait(msgbox(message, 'Cleaning params changed'));

%% Set up parameters for AFQ_AddNewFiberGroup

% List of ROI names for the different SLF branches
% ROIs go frontal to parietal

roi1Names = {...
    'SFgL.nii.gz' 'SFgR.nii.gz'...
    'MFgL.nii.gz' 'MFgR.nii.gz'...
    'PrgL.nii.gz' 'PrgR.nii.gz'
    };

roi2Names = {...
    'PaL.nii.gz' 'PaR.nii.gz'...
    'PaL.nii.gz' 'PaR.nii.gz'...
    'PaL.nii.gz' 'PaR.nii.gz'
    };

% Temporal 'NOT' ROI to exclude frontotemporal fibers of the arcuate
% fascilulus
roi3Names = {...
    'SLFt_roi2_L.nii.gz' 'SLFt_roi2_R.nii.gz'...
    'SLFt_roi2_L.nii.gz' 'SLFt_roi2_R.nii.gz'...
    'SLFt_roi2_L.nii.gz' 'SLFt_roi2_R.nii.gz'
    };

% Specify path to SLF123 ROIs
tdir = uigetdir('','Choose the directory where the SLF ROIs are stored');
for i = 1:length(roi1Names)
    roi1fullNames{i} = fullfile(tdir,roi1Names{i});
    roi2fullNames{i} = fullfile(tdir,roi2Names{i});
    roi3fullNames{i} = fullfile(tdir,roi3Names{i});
end

% List of the SLF branches
fgNames = {...
    'LeftSLF1.mat' 'RightSLF1.mat'...
    'LeftSLF2.mat' 'RightSLF2.mat'...
    'LeftSLF3.mat' 'RightSLF3.mat'
    };


%% Segment SLF into its different projections and add to afq structure
segNumber=length(afq.fgnames);

for i = 1:length(fgNames)
    afq = SLF123_AFQ_AddNewFiberGroup3ROI(afq,fgNames{i},roi1fullNames{i},roi2fullNames{i},roi3fullNames{i},'not',1,1,0,[],overwriteFiles,segNumber+i);
end

% Change afq cleaning params to original values:
afq.params.maxDist = orig_dist;
afq.params.maxLen = orig_len;

save(afq_fileName,'afq');
end

