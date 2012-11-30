function afq = AFQ_AddNewFiberGroup(afq,fgName,roi1Name,roi2Name,cleanFibers,computeVals)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Add new fiber groups from any segmentation proceedure to an afq structure
%
% afq = AFQ_AddNewFiberGroup(afq, fgNaem, roi1Name, roi2Name, [cleanFibers = true],[computeVals = true])
%
% By default AFQ_run will segment the 20 fiber groups defined in the Mori
% white matter atlas and save all the relevant information in the afq
% structure. Aditional groups defined manually or with other segmentation
% algorithms can be added to the afq structure with AFQ_AddNewFiberGroup.
% By default this function will take in (1) an afq strucutre that has already
% been passed through AFQ_run (see AFQ_Create and AFQ_run), (2) the name of
% the new fiber group which should be in each subject's fibers directory,
% and (3) the name of two waypoint ROIs that the fiber group passes through
% which should be in each subjects ROIs directory. The necesary fields in
% the afq structure will be populated, TractProfiles will be calculated and
% added to the structure (See AFQ_CreateTractProfile and
% AFQ_ComputeTractProperties) and the new afq structure will be returned.
% Feel free to take an already segmented fiber group and associate it with
% 2 new ROIs in order to define the edges of the profile differently.
%
% Inputs:
% afq        - An afq structure output from AFQ_run (see also AFQ_Create)
% fgName     - A string containing the name of the fiber group. The fiber
%              group should be in each subject's fibers directory (with the
%              same name). This can be either a .mat file or a .pdb file.
%              If the fgName file does not exist then the fiber group will
%              be defined by intersecting the subject's wholeBrainFG with
%              the two ROIs.
% roi1Name   - A string containing the name of the first waypoint ROI.
%              The ROI should be in each subject's ROIs directory
% roi2Name   - A string containing the name of the second ROI
% cleanFibers- Locigal defining whether the fibers should be cleaned. See
%              AFQ_removeFiberOutliers
% computeVals- Logical defining whether TractProfiles should be calculated
%              for the new fiber group for all of the diffusion (and other)
%              parameters that are  already contained in the afq structure
%
% Copyright Jason D. Yeatman November 2012

%% Argument checking
if ~isafq(afq)
    error('Please enter an afq structure')
end
if ~exist('fgName','var') || isempty(fgName) || ~ischar(fgName)
    error('Please enter the name of the fiber group')
elseif ~strcmp(fgName(end-3:end),'.mat') || ~strcmp(fgName(end-3:end),'.pdb')
    % add file extension
    fgName = [fgName '.mat'];
end
if ~exist('roi1Name','var') || isempty(roi1Name) || ~ischar(roi1Name)
    error('Please enter the name of the first ROI')
elseif ~strcmp(roi1Name(end-3:end),'.mat')
    roi1Name = [roi1Name '.mat'];
end
if ~exist('roi2Name','var') || isempty(roi2Name) || ~ischar(roi2Name)
    error('Please enter the name of the second ROI')
elseif ~strcmp(roi2Name(end-3:end),'.mat')
    roi2Name = [roi2Name '.mat'];
end
if ~exist('cleanFibers','var') || isempty(cleanFibers)
    cleanFibers = true;
end
if ~exist('computeVals','var') || isempty(computeVals)
    computeVals = true;
end

%% Add the new fiber groups and ROIs to the afq structure
afq = AFQ_set(afq,'new fiber group', fgName);
afq = AFQ_set(afq, 'new roi', roi1Name, roi2Name);
% Get the fiber group number. This will be equal to the number of fiber
% groups since it is the last one to be added
fgNumber = AFQ_get(afq,'numfg');

%% Segment the fiber groups if they don't exist
for ii = 1:AFQ_get(afq,'numsubs')
    if ~exist(AFQ_get(afq,[prefix(fgName) 'path'],ii),'file')
        % Load the wholebrain fiber group
        wholebrainFG = AFQ_get(afq,'wholebrain fg', ii);
        % Load the defining ROIs
        [roi1, roi2] = AFQ_LoadROIs(fgNumber,afq.sub_dir{ii}, afq);
        % Intersect the wholebrain fibers with each ROI
        fg_classified = dtiIntersectFibersWithRoi([],'and', 2, roi1, wholebrainFG);
        fg_classified = dtiIntersectFibersWithRoi([],'and', 2, roi2, fg_classified);
        % Set the name
        fg_classified.name = prefix(fgName);
        % Save it
        dtiWriteFiberGroup(fg_classified,AFQ_get(afq,[prefix(fgName) 'path'],ii));
    end
end

%% Clean the fibers if desired
if cleanFibers == 1
    for ii 1::AFQ_get(afq,'numsubs')
        % Load the fibers
        if ~exist('fg_classified','var')
            fg_classified = dtiLoadFiberGroup(AFQ_get(afq,[prefix(fgName) 'path'],ii));
        end
        
        % only clean if there are enough fibers for it to be worthwhile
        if  length(fg_classified.fibers) > 20
            % clean clipped fiber group if computations are to be done
            % on clipped group
            if afq.params.cleanClippedFibers == 1;
                % load ROIs
                [roi1, roi2] = AFQ_LoadROIs(jj,afq.sub_dirs{ii},afq);
                % clip fiber group
                fg_clip      = dtiClipFiberGroupToROIs(fg_classified,roi1,roi2);
                if length(fg_clip.fibers) > 20
                    % clean clipped fiber group
                    [~, keep] = AFQ_removeFiberOutliers(fg_classified,afq.params.maxDist,afq.params.maxLen,afq.params.numberOfNodes,'mean',0);
                    % remove fibers from unclipped group that do no
                    % survive the clipping
                    fg_classified = fg_classified.fibers(keep);
                end
            else
                % clean un-clipped fiber group
                fg_classified = AFQ_removeFiberOutliers(fg_classified,afq.params.maxDist,afq.params.maxLen,afq.params.numberOfNodes,'mean',0,afq.params.cleanIter);
            end
        end
        
        % Define the full path to the new cleaned fiber group
        fgpath = fullfile(afq.sub_dirs{ii},'fibers',[prefix(fgName) '_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']);
        % Save the fiber group
        dtiWriteFiber(fg, fgpath);
        % And add them to the afq structure
        afq.files.fibers.([prefix(fgName) 'clean']){ii} = fgpath;
    end
end

%% Compute tract profiles