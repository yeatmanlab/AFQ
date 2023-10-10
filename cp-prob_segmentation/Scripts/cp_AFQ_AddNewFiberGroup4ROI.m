function afq = cp_AFQ_AddNewFiberGroup4ROI(afq,fgName,roi1Name,roi2Name,roi3Name,roi3type,roi4Name,roi4type, cleanFibers,computeVals,showFibers,segFgName,overwrite,fgNumber)
% This function is based on cp_AFQ_AddNewFiberGroup3ROI:

% September 2021 - Sivan Jossinger
% Added roi4Name to easily add fibergroups with an additional AND roi, and NOT roi, before cleaning
% stage

% 17 August 2015 - Maya Yablonski
% Added roi3Name to easily add fibergroups with a NOT roi, before cleaning
% stage

% New inputs:

% roi3Name- name of additional ROI to intersect the fiber with
% roi3type- defines the kind of intersection with the third ROI e.g. 'and','not'

% roi4Name- name of additional ROI to intersect the fiber with
% roi4type- defines the kind of intersection with the third ROI e.g. 'and','not'


% afq = AFQ_AddNewFiberGroup(afq, fgName, roi1Name, roi2Name, roi3Name,roi3type, roi3Name,roi3type, [cleanFibers = true], ...
%          [computeVals = true], [showFibers = false], [segFgName = 'WholeBrainFG.mat'] ...
%          [overwrite = false])
%

%
% Inputs:
% afq        - An afq structure output from AFQ_run (see also AFQ_Create)
% fgName     - A string containing the name of the fiber group. The fiber
%              group should be in each subject's fibers directory (with the
%              same name). This can be either a .mat file or a .pdb file.
%              If the fgName file does not exist then the fiber group will
%              be defined by intersecting the subject's wholeBrainFG with
%              the two ROIs.
% roi1Name   - A string containing the name of the first waypoint ROI. This
%              can either be a .mat file that is in each subject's ROIs
%              directory or it can be a nifti image containing the ROI
%              defined in MNI space. If it is a nifti image than that image
%              will be warped into each subject's native space and saved as
%              a .mat ROI. This allows new segmentations to be added to AFQ
%              simply by defining the ROIs on the MNI template
% roi2Name   - A string containing the name of the second ROI (either .mat
%              files that are in each subject's directory or a single nifti
%              image in MNI space)
% roi3Name   - A string containing the name of the third ROI (either .mat
%              files that are in each subject's directory or a single nifti
%              image in MNI space). 
% roi3Type   - The type of intersectiuon with the additional ROI, e.g.,
%               'and, 'not' etc.  
% roi4Name   - A string containing the name of the third ROI (either .mat
%              files that are in each subject's directory or a single nifti
%              image in MNI space). 
% roi4Type   - The type of intersectiuon with the additional ROI, e.g.,
%               'and, 'not' etc.  
% cleanFibers- Logical defining whether the fibers should be cleaned. See
%              AFQ_removeFiberOutliers
% computeVals- Logical defining whether TractProfiles should be calculated
%              for the new fiber group for all of the diffusion (and other)
%              parameters that are  already contained in the afq structure
% showFibers - Logical indicating whether to render each segmented and
%              cleaned fiber group.
% segFgName  - There are some case where you may want to perform the
%              segmentation on a fiber group other than the wholebrain
%              group (eg. for the callosum). If you would like to use
%              another fiber group you can supply it's name here otherwise
%              WholeBrainFG.mat will be used
% overwrite  - Whether or not to overwrite previously computed files
%
%
% Copyright Jason D. Yeatman November 2012

%% Argument checking
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end
if ~isafq(afq)
    error('Please enter an afq structure')
end
if ~exist('fgName','var') || isempty(fgName) || ~ischar(fgName)
    error('Please enter the name of the fiber group')
elseif ~(strcmp(fgName(end-3:end),'.mat') || strcmp(fgName(end-3:end),'.pdb'))
    % add file extension
    fgName = [fgName '.mat'];
end
% Check that the roi names were defined correctly
if ~exist('roi1Name','var') || isempty(roi1Name) || ~ischar(roi1Name)
    error('Please enter the name of the first ROI')
elseif ~isempty(strfind(roi1Name,'.nii'))
    % If the roi names were nifti images then we are assuming they are
    % defined in MNI space and should be transformed to each individual's
    % brain
    if exist(roi1Name,'file') && exist(roi2Name,'file') && exist(roi3Name,'file') && exist(roi4Name,'file')
        xformRois = 1;
        fprintf('\nAssuming the ROI is defined in MNI space. It will be xformed to each individuals brain\n')
        % Assign roi names to template roi variables
        Troi1 = roi1Name; Troi2 = roi2Name;
        % MY - same for third ROI
        Troi3 = roi3Name;
        % SJ - same for third ROI
        Troi4 = roi4Name;    
        % Get the names of the rois without paths and .nii suffix
        [~,tmp] = fileparts(roi1Name); roi1Name = [prefix(tmp) '.mat'];
        [~,tmp] = fileparts(roi2Name); roi2Name = [prefix(tmp) '.mat'];
        [~,tmp] = fileparts(roi3Name); roi3Name = [prefix(tmp) '.mat'];
        [~,tmp] = fileparts(roi4Name); roi4Name = [prefix(tmp) '.mat'];
    else
        error('Could not find template ROI file')
    end
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
if ~exist('showFibers','var') || isempty(showFibers)
    showFibers = false;
end
if ~exist('segFgName','var') || isempty(segFgName)
    segFgName = 'WholeBrainFG.mat';
end
% If the ROIs were not defined as a nifti image in MNI space than they do
% not need to be transformed
if ~exist('xformRois','var') || isempty(xformRois)
    xformRois = 0;
end
% Check which subjects should be run
runsubs = AFQ_get(afq,'run subjects');

%% Add the new fiber groups and ROIs to the afq structure
if notDefined('fgNumber') || fgNumber > AFQ_get(afq,'numfg');
    afq = AFQ_set(afq,'new fiber group', fgName);
    afq = AFQ_set(afq, 'new roi', roi1Name, roi2Name);
    % Get the fiber group number. This will be equal to the number of fiber
    % groups since it is the last one to be added
    fgNumber = AFQ_get(afq,'numfg');
end

%% Make individual ROIs from a template ROI if a template was passed in
if xformRois == 1
    % Path to the templates directory
    tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');  
    % template = fullfile(tdir,'MNI_JHU_T2.nii.gz');
    template = AFQ_get(afq, 'template'); % Jason fix 19/07/2016 %
        
    for ii = runsubs
        % Get the subject's dt6 file
        dtpath = AFQ_get(afq,'dt6path',ii); dt = dtiLoadDt6(dtpath);
        % Subject's directory
        sdir = fileparts(dtpath);
        % Check if the ROIs already exist
        if      ~exist(fullfile(sdir,'ROIs',roi1Name),'file') || ...
                ~exist(fullfile(sdir,'ROIs',roi2Name),'file') ||...
                ~exist(fullfile(sdir,'ROIs',roi3Name),'file') ||...
                ~exist(fullfile(sdir,'ROIs',roi4Name),'file') ||...
            overwrite == 1    
            % prevline MY
            % Check if there is a precomputed spatial normalization.
            % Otherwise compute spatial normalization
            sn = AFQ_get(afq,'spatial normalization',ii);
            invDef = AFQ_get(afq,'inverse deformation',ii);
            if isempty(sn) || isempty(invDef)
                %[sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(dt.b0, dt.xformToAcpc, template);
                % Rescale image valueds to get better gray/white/CSF contrast
                alignIm = mrAnatHistogramClip(double(dt.b0),0.3,0.99);
                [sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(alignIm, dt.xformToAcpc, template);
            end          
            % Load up template ROIs in MNI space and transform them to the subjects
            % native space.
            [~, ~, roi1]=dtiCreateRoiFromMniNifti(dt.dataFile, Troi1, invDef, 0);
            [~, ~, roi2]=dtiCreateRoiFromMniNifti(dt.dataFile, Troi2, invDef, 0);
            [~, ~, roi3]=dtiCreateRoiFromMniNifti(dt.dataFile, Troi3, invDef, 0);
            [~, ~, roi4]=dtiCreateRoiFromMniNifti(dt.dataFile, Troi4, invDef, 0);

            % same for third ROI
            % Save the ROIs as .mat files. Don't overwrite if roi already exists in subject's directory
         
            if ~exist(fullfile(sdir,'ROIs',roi1Name),'file')
            dtiWriteRoi(roi1,fullfile(sdir,'ROIs',roi1Name));
            end
            if ~exist(fullfile(sdir,'ROIs',roi2Name),'file')
            dtiWriteRoi(roi2,fullfile(sdir,'ROIs',roi2Name));
            end
            if ~exist(fullfile(sdir,'ROIs',roi3Name),'file')
            dtiWriteRoi(roi3,fullfile(sdir,'ROIs',roi3Name));
            end
            if ~exist(fullfile(sdir,'ROIs',roi4Name),'file')
            dtiWriteRoi(roi4,fullfile(sdir,'ROIs',roi4Name));
            end
        else
            fprintf('\n %s, %s, %s and %s exist for subject %d',roi1Name,roi2Name,roi3Name,roi4Name,ii)
        end
    end
end

%% Segment the fiber groups if they don't exist
for ii = runsubs
    
    % Define the current subject to process
    afq = AFQ_set(afq,'current subject',ii);
    
    if ~exist(AFQ_get(afq,[prefix(fgName) 'path'],ii),'file') || overwrite == 1
        % Load the wholebrain fiber group as default
        % or use another fiber group if desired (eg. callosum)
        segFgPath = fullfile(afq.sub_dirs{ii}, 'fibers', segFgName);
        if exist(segFgPath,'file')
            fprintf('\nPerforming segmentation on %s',segFgPath)
            wholebrainFG = dtiLoadFiberGroup(segFgPath);
        else
            error('\nCould not find %s',segFgPath);
        end

        % Load the defining ROIs
        [roi1, roi2] = AFQ_LoadROIs(fgNumber,afq.sub_dirs{ii}, afq);
        
        % use AFQ_get to get current subject directory sdir
        dtpath = AFQ_get(afq,'dt6path',ii);
        sdir = fileparts(dtpath);          
        roi3 = dtiReadRoi(fullfile(sdir,'ROIs',roi3Name));
        roi4 = dtiReadRoi(fullfile(sdir,'ROIs',roi4Name));
        % Intersect the wholebrain fibers with each ROI
        fg_classified = dtiIntersectFibersWithRoi([],'and', 2, roi1, wholebrainFG);
        fg_classified = dtiIntersectFibersWithRoi([],'and', 2, roi2, fg_classified);
        % To remove midline fibers, intersect with roi3 using the default
        % minDist (empty array, dtiIntersectFibersWithRoi will use the
        % default 0.87). For MCP segmentation, this was found to give
        % better results than minDist 2. -M.Y.
        % Consider using minDist 2 if using this function
        % for other fibers in the future. 
        %fg_classified = dtiIntersectFibersWithRoi([],roi3type, 2, roi3, fg_classified);
        fg_classified = dtiIntersectFibersWithRoi([],roi3type, 2, roi3, fg_classified);
        fg_classified = dtiIntersectFibersWithRoi([],roi4type, 2, roi4, fg_classified);

        % MY August 2015 - add intersection with third ROI
        % Flip the fibers so they all pass through roi1 before roi2
        fg_classified = AFQ_ReorientFibers(fg_classified,roi1,roi2);
        % Set the name
        fg_classified.name = prefix(fgName);
        % Save it
        dtiWriteFiberGroup(fg_classified,AFQ_get(afq,[prefix(fgName) 'path'],ii));
        % Clear variables
        clear fg_classified wholebrainFG roi1 roi2 roi3 roi4
    end
end

%% Clean the fibers if desired
for ii = runsubs
    
    % Define the current subject to process
    afq = AFQ_set(afq,'current subject',ii);
    % Define the full path to the new cleaned fiber group
    fgclean_path = fullfile(afq.sub_dirs{ii},'fibers',[prefix(fgName) '_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']);
    
    % Only clean if desired and if the cleaned fibers do not already exist
    if cleanFibers == 1 && (~exist(fgclean_path,'file') || overwrite == 1)
        
        % Get path to fibers
        fg_classified_path = AFQ_get(afq,[prefix(fgName) 'path'],ii);
        fprintf('\nCleaning %s',fg_classified_path);
        % Load the fibers
        fg_classified = dtiLoadFiberGroup(fg_classified_path);
        
        % Render the segmented fibers if desired
        if showFibers == 1
            fprintf('\n Rendering %s in red\n',AFQ_get(afq,[prefix(fgName) 'path'],ii));
            AFQ_RenderFibers(fg_classified, 'numfibers',70,'color',[1 0 0])
        end
        
        % only clean if there are enough fibers for it to be worthwhile
        if  length(fg_classified.fibers) > 20
            % clean clipped fiber group if computations are to be done
            % on clipped group
            if afq.params.cleanClippedFibers == 1;
                % load ROIs
                [roi1, roi2] = AFQ_LoadROIs(jj,afq.sub_dirs{ii},afq);
                % clip fiber group
                % fg_clip      = dtiClipFiberGroupToROIs(fg_classified,roi1,roi2);
                fg_clip      = cp_dtiClipFiberGroupToROIs(fg_classified,roi1,roi2);
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
        
        % Save the fiber group
        dtiWriteFiberGroup(fg_classified, fgclean_path);
        % And add them to the afq structure
        afq.files.fibers.([prefix(fgName) '_clean']){ii} = fgclean_path;
        
        % Render the cleaned fibers if desired
        if showFibers == 1
            fprintf('\n Rendering %s in blue\n',AFQ_get(afq,[prefix(fgName) 'path'],ii));
            AFQ_RenderFibers(fg_classified, 'numfibers',70,'color',[0 0 1])
        end  
        % clear the variables
        clear fg_classified fg_classified_path fgclean_path
        
    elseif exist(fgclean_path,'file')
        fprintf('\n%s already exists',fgclean_path);
        % Add the path to the cleaned fiber group to the afq structure
        afq.files.fibers.([prefix(fgName) '_clean']){ii} = fgclean_path;
    end
end

%% Compute tract profiles
for ii = runsubs
    
    % Define the current subject to process
    afq = AFQ_set(afq,'current subject',ii);
    
    if computeVals == 1
        fprintf('\nComputing Tract Profiles for subject %s',afq.sub_dirs{ii});
        % Load the subject's dt6
        dt = dtiLoadDt6(AFQ_get(afq,'dt6path',ii));
        % Load the propper fiber group. If they have been cleaned it will
        % load the cleaned version
        fg_classified_name = AFQ_get(afq,[prefix(fgName) 'path'],ii);
        fprintf('\nFiber group: %s',fg_classified_name);
        fg_classified = dtiReadFibers(fg_classified_name);
        
        % Determine how much to weight each fiber's contribution to the
        % measurement at the tract core. Higher values mean steaper falloff
        fWeight = AFQ_get(afq,'fiber weighting');
        % By default Tract Profiles of diffusion properties will always be
        % calculated
        %[fa, md, rd, ad, cl, volume, TractProfile] = AFQ_ComputeTractProperties(fg_classified, dt, afq.params.numberOfNodes, afq.params.clip2rois, afq.sub_dirs{ii}, fWeight, afq);
        [fa, md, rd, ad, cl, volume, TractProfile] = cp_AFQ_ComputeTractProperties(fg_classified, dt, afq.params.numberOfNodes, afq.params.clip2rois, afq.sub_dirs{ii}, fWeight, afq);
        
        % Parameterize the shape of each fiber group with calculations of
        % curvature and torsion at each point and add it to the tract
        % profile
        [curv, tors, TractProfile] = AFQ_ParamaterizeTractShape(fg_classified, TractProfile);
        % Fill these variables with nans if they come out empty. This will
        % happen if the tract profile is empty for this subject
        if isempty(curv) || isempty(tors)
            curv = nan(size(fa)); tors = nan(size(fa));
        end
        
        % Calculate the volume of each Tract Profile
        TractProfile = AFQ_TractProfileVolume(TractProfile);
        
        % Add values to the afq structure
        afq = AFQ_set(afq,'vals','subnum',ii,'fgnum',fgNumber, 'fa',fa, ...
            'md',md,'rd',rd,'ad',ad,'cl',cl,'volume',volume,'curvature',curv,'torsion',tors);
        
        % Add Tract Profiles to the afq structure
        afq = AFQ_set(afq,'tract profile','subnum',ii,'fgnum',fgNumber,TractProfile);
        
        % If any other images were supplied calculate a Tract Profile for that
        % parameter
        numimages = AFQ_get(afq, 'numimages');
        if numimages > 0;
            for jj = 1:numimages
                % Read the image file
                image = readFileNifti(afq.files.images(jj).path{ii});
                % Compute a Tract Profile for that image
                imagevals = cp_AFQ_ComputeTractProperties(fg_classified, image, afq.params.numberOfNodes, afq.params.clip2rois, afq.sub_dirs{ii}, fWeight, afq);
                % Add values to the afq structure
                afq = AFQ_set(afq,'vals','subnum',ii,'fgnum',fgNumber, afq.files.images(jj).name, imagevals);
                clear imagevals
            end
        end
    else
        fprintf('\nTract Profiles already computed for subject %s',afq.sub_dirs{ii});
    end
    
    % Save each iteration of afq run if an output directory was defined
    if ~isempty(AFQ_get(afq,'outdir')) && exist(AFQ_get(afq,'outdir'),'dir')
        if ~isempty(AFQ_get(afq,'outname'))
            outname = fullfile(AFQ_get(afq,'outdir'),AFQ_get(afq,'outname'));
        else
            outname = fullfile(AFQ_get(afq,'outdir'),['afq_' date]);
        end
        save(outname,'afq');
    end
end

%% Recompute the norms with the new fiber group
if AFQ_get(afq,'computenorms')
    [norms, patient_data, control_data, afq] = AFQ_ComputeNorms(afq);
end



