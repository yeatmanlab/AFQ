function afq = AFQ_SegmentCallosum(afq,overwriteFiles,sge)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Segment the callosal fibers into 7 segments
%
% afq = AFQ_SegmentCallosum(afq,overwriteFiles)
%
%
%
% Copyright Jason D. Yeatman April 2012

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
%% Create callosum ROIs and fiber groups
% First creat a callosal fiber group from the wholebrain fiber group for
% each subject
for ii = 1:AFQ_get(afq,'numsubs')
    % Set up paths to the fiber group and ROI
    fgPath = fullfile(afq.sub_dirs{ii},'fibers','callosumFG.mat');
    roiPath = fullfile(afq.sub_dirs{ii},'ROIs','callosum_rough.mat');
    if ~exist(fgPath,'file') || ~exist(roiPath,'file') || overwriteFiles(1)==1
        % Load Dt6
        dt = dtiLoadDt6(AFQ_get(afq,'dt6path',ii));
        % Create an ROI of the corpus callosum
        %ccRoi = dtiNewRoi('callosum','r',dtiFindCallosum(dt.dt6,dt.b0,dt.xformToAcpc,.25,[],1));
        [~,~,ccRoi]=dtiCreateRoiFromMniNifti(dt.dataFile, fullfile(AFQ_directories,'templates','callosum','callosum_rough.nii.gz'));
        % Load wholebrain fiber group
        wholebrainFg = AFQ_get(afq,'wholebrain fg', ii);
        % Create callosum fiber group
        ccFg = dtiIntersectFibersWithRoi([],'and',2,ccRoi,wholebrainFg);
        ccFg.name = 'callosumFG';
        % Save callosum ROI and fiber group
        fprintf('\nSaving %s',fgPath)
        dtiWriteFiberGroup(ccFg,fgPath);
        dtiWriteRoi(ccRoi,roiPath);
    end
end

%% Set up parameters for AFQ_AddNewFiberGroup

% List of ROI names for different callosal segments
roi1Names = {'L_Occipital.nii.gz' 'L_PostParietal.nii.gz' ...
    'L_SupParietal.nii.gz' 'L_Motor.nii.gz' 'L_SupFrontal.nii.gz' ...
    'L_AntFrontal.nii.gz' 'L_Orbital.nii.gz' 'L_Temporal.nii.gz'};
roi2Names = {'R_Occipital.nii.gz' 'R_PostParietal.nii.gz' ...
    'R_SupParietal.nii.gz' 'R_Motor.nii.gz' 'R_SupFrontal.nii.gz' ...
    'R_AntFrontal.nii.gz' 'R_Orbital.nii.gz' 'R_Temporal.nii.gz'};
% Add in the path to the templates folder to the roi paths
tdir = fullfile(AFQ_directories,'templates','callosum2');
for ii = 1:length(roi1Names)
    roi1Names{ii} = fullfile(tdir,roi1Names{ii});
    roi2Names{ii} = fullfile(tdir,roi2Names{ii});
end

% List of the names of the callosal segments
fgNames = {'CC_Occipital.mat' 'CC_Post_Parietal.mat' ...
    'CC_Sup_Parietal.mat' 'CC_Motor.mat' 'CC_Sup_Frontal.mat' ...
    'CC_Ant_Frontal.mat' 'CC_Orb_Frontal.mat' 'CC_Temporal.mat'};
% The name of the fiber group file to segment
segFgName = 'callosumFG.mat';

%% Segment callosum into it's different projections and add to afq struct

% Run using sun grid engine if desired
if sge == 1
    for ii = 1:length(fgNames)
        afq = AFQ_AddNewFiberGroup_sge(afq,fgNames{ii},roi1Names{ii},roi2Names{ii},1,1,0,segFgName,overwriteFiles(2));
    end
else
    for ii = 1:length(fgNames)
        afq = AFQ_AddNewFiberGroup(afq,fgNames{ii},roi1Names{ii},roi2Names{ii},1,1,0,segFgName,overwriteFiles(2),20+ii);
    end
end

%% Render montage of callosal fiber groups
fgNames = {'CC_Occipital' 'CC_Post_Parietal' ...
    'CC_Sup_Parietal' 'CC_Motor' 'CC_Sup_Frontal' ...
    'CC_Ant_Frontal' 'CC_Orb_Frontal' 'CC_Temporal'};
AFQ_MakeFiberGroupMontage(afq, fgNames)

return