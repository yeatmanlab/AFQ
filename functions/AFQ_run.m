function [afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq)
% Run AFQ analysis on a set of subjects to generate Tract Profiles of white
% matter properties.
%
% [afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs,
% sub_group, [afq])
%
% AFQ_run is the main function to run the AFQ analysis pipeline.  Each AFQ
% function is an independent module that can be run on its own.  However
% when AFQ_run is used to analyze data all the proceeding analyses are
% organized into the afq data structure. The AFQ analysis pipeline is
% described in Yeatman J.D., Dougherty R.F., Myall N.J., Wandell B.A.,
% Feldman H.M. (2012). Tract Profiles of White Matter Properties:
% Automating Fiber-Tract Quantification. PLoS One.
%
% Input arguments:
%  sub_dirs  = 1 x N cell array where N is the number of subjects in the
%              study. Each cell should contain the full path to a subjects
%              data directory where there dt6.mat file is.
%
%  sub_group = Binary vector defining each subject's group. 0 for control
%              and 1 for patient.
%
%  afq       = This is a structure that sets up all the parameters for the
%              analysis.  If it is blank AFQ_run will use the default
%              parameters.  See AFQ_Create.
%
% Outputs: 
% afq          = afq structure containing all the results
%
% patient_data = A 1X20 structured array of tract diffusion profiles where
%                data for each tract is in a cell of the structure (eg.
%                patient_data(1) is data for the left thalamic radiation).
%                Each diffusion properties is stored as a different field
%                (eg. patient_data(1).FA is a matrix of FA profiles for the
%                left thalamic radiation). Within the data matrix each
%                subject is a row and each location is a column.  This
%                output variable contains all the data for the patients
%                defined by sub_group ==1.
%
% control_data = The same structure as for patient_data but this contains
%                data for the control subjects defined by sub_group==0.
%
% norms        = Means and standard deviations for each tract diffusion
%                profile calculated based on the control_data.
%
% abn          = A 1 x N vector where N is the number of patients.
%                Each patient that is abnormal on at least one tract is
%                marked with a 1 and each subject that is normal on every
%                tract is marked with a 0. The criteria for abnormal is
%                defined in afq.params.cutoff.  See AFQ create
%
% abnTracts    = An M by N matrix where M is the number of subjects and N
%                is the number of tracts. Each row is a subject and each 
%                column is a tract.  1 means that tract was abnormal for
%                that subject and 0 means it was normal.
%
%  Web resources
%    http://white.stanford.edu/newlm/index.php/AFQ
%
%  Example:
%   
%   % Get the path to the AFQ directories
%   [AFQbase AFQdata] = AFQ_directories;
%   % Create a cell array where each cell is the path to a data directory
%   sub_dirs = {[AFQdata '/patient_01/dti30'], [AFQdata '/patient_02/dti30']...
%   [AFQdata '/patient_03/dti30'], [AFQdata '/control_01/dti30']...
%   [AFQdata '/control_02/dti30'], [AFQdata '/control_03/dti30']};
%   % Create a vector of 0s and 1s defining who is a patient and a control
%   sub_group = [1, 1, 1, 0, 0, 0]; 
%   % Run AFQ in test mode to save time. No inputs are needed to run AFQ 
%   % with the default settings. AFQ_Create builds the afq structure. This
%   % will also be done automatically by AFQ_run if the user does not wish 
%   % to modify any parameters
%   afq = AFQ_Create('run_mode','test', 'sub_dirs', sub_dirs, 'sub_group', sub_group); 
%   [afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq)
%
% Copyright Stanford Vista Team, 2011. Written by Jason D. Yeatman,
% Brian A. Wandell and Robert F. Dougherty

%% Check Inputs
if notDefined('sub_dirs') && exist('afq','var') && ~isempty(afq)
    sub_dirs = AFQ_get(afq,'sub_dirs');
elseif notDefined('sub_dirs')
    error('No subject directories');
end
if ~iscell(sub_dirs), sub_dirs = cellstr(sub_dirs); end
if notDefined('sub_group') && exist('afq','var') && ~isempty(afq)
    sub_group = AFQ_get(afq,'sub_group');
elseif notDefined('sub_group')
    error('Must define subject group');
end
if length(sub_group) ~= size(sub_dirs,1) && length(sub_group) ~= size(sub_dirs,2)
    error('Mis-match between subject group description and subject data directories');
end
if ~exist('afq','var') || isempty(afq)
    %if no parameters are defined use the defualts
    afq = AFQ_Create('sub_dirs',sub_dirs,'sub_group',sub_group); 
end
if isempty(afq.sub_group)
    afq = AFQ_set(afq,'sub_group',sub_group);
end
if isempty(afq.sub_dirs)
    afq = AFQ_set(afq,'sub_dirs',sub_dirs);
end
% Check which subjects should be run
runsubs = AFQ_get(afq,'run subjects');
% Define the name of the segmented fiber group
segName = AFQ_get(afq,'segfilename');

% If ANTS is installed on the system then precompute spatial normalization
% with ANTS and save to the afq structure
if AFQ_get(afq, 'use ANTS')
    afq = AFQ_ComputeSpatialNormalization(afq);
end

%%  Loop over every subject
for ii = runsubs
    % Define the current subject to process
    afq = AFQ_set(afq,'current subject',ii);
    
    %% Preprocess Data
    
    if ~exist(fullfile(sub_dirs{ii},'dt6.mat'),'file')
        error('AFQ preprocessing has not yet been implemented')
    end
    % Load Dt6 File
    dtFile = fullfile(sub_dirs{ii},'dt6.mat');
    dt     = dtiLoadDt6(dtFile);
    
    % If ANTS was used to compute a spatial normalization then load it for
    % this subject
    antsInvWarp = AFQ_get(afq,'ants inverse warp',ii);
    
    %% Perform Whole Brain Streamlines Tractography
    
    %Check if there is a fibers directory, otherwise make one.
    fibDir = fullfile(sub_dirs{ii},'fibers');
    if ~exist(fibDir,'dir')
        mkdir(fibDir);
    end 
    % Check if wholebrain tractography should be done
    if AFQ_get(afq, 'do tracking',ii) == 1
        % Perform whole-brain tractography
        fprintf('\nPerforming whole-brain tractograpy for subject %s\n',sub_dirs{ii});
        fg = AFQ_WholebrainTractography(dt, afq.params.run_mode, afq);
        % Save fiber group to fibers directory
        dtiWriteFiberGroup(fg,fullfile(fibDir,'WholeBrainFG.mat'));
        % Set the path to the fibers in the afq structure
        afq = AFQ_set(afq,'wholebrain fg path','subnum',ii,fullfile(fibDir,'WholeBrainFG.mat'));
        % Wholebrain fiber group is already in memory and does not need to
        % be loaded
        loadWholebrain = 0;
    else
        fprintf('\nWhole-brain tractography was already done for subject %s',sub_dirs{ii});
        % Wholebrain fiber group needs to be loaded
        loadWholebrain = 1;
    end
    
    %% Segment 20 Fiber Groups
    
    % Check if fiber group segmentation was already done
    if AFQ_get(afq, 'do segmentation',ii) == 1
        % Load the wholebrain fiber group if necessary
        if loadWholebrain == 1
            fg = AFQ_get(afq,'wholebrain fiber group',ii);
        end
        % Segment fiber group
        fg_classified = AFQ_SegmentFiberGroups(dtFile, fg, [], [],[], antsInvWarp);
        % Save segmented fiber group
        dtiWriteFiberGroup(fg_classified, fullfile(fibDir,segName));
        % If the full trajectory of each fiber group will be analyzed (eg.
        % from cortical start to endpoint) then all fibers that terminate
        % before cortex will be removed and each fiber within a group will
        % be flipped so that startpoints and endpoints are consistent
        % within the group
        if AFQ_get(afq,'clip2rois') == 0
            fg_classified = AFQ_DefineFgEndpoints(fg_classified, [], [], dt);
            dtiWriteFiberGroup(fg_classified, fullfile(fibDir,segName));
        end
        % Set the path to the fibers in the afq structure
        afq = AFQ_set(afq, 'segmented fg path', 'subnum', ii, fullfile(fibDir,segName));
        % Segemented fiber group is already in memory and does not need to
        % be loaded
        loadSegmentation = 0;
        clear fg
    else
        fprintf('\nFiber tract segmentation was already done for subject %s',sub_dirs{ii});
        % Segmentation needs to be loaded
        loadSegmentation = 1;
    end
    clear dtFile fgFile
    
    %% Remove fiber outliers from each fiber tract so it is a compact bundle    
    
    % Run AFQ cleaning proceedure if it is called for in
    % afq.params.cleanFibers and if has not yet been done
    if afq.params.cleanFibers == 1 && AFQ_get(afq, 'do cleaning', ii) == 1
        % Load segmented fiber group if necessary
        if loadSegmentation == 1
            fg_classified = dtiLoadFiberGroup(fullfile(fibDir, segName));
        end
        % Convert fiber groups into an array if they are not already 
        fg_clean = fg2Array(fg_classified); 
        
        % Remove all fibers that are too long and too far from the core of
        % the group.  This algorithm will constrain the fiber group to
        % something that can be reasonable represented as a 3d gaussian
        for jj = 1:20
            % only clean if there are enough fibers for it to be worthwhile
            if  length(fg_clean(jj).fibers) > 20
                % clean clipped fiber group if computations are to be done
                % on clipped group
                if afq.params.cleanClippedFibers == 1;
                    % load ROIs
                    [roi1 roi2] = AFQ_LoadROIs(jj,sub_dirs{ii});
                    % clip fiber group
                    fg_clip     = dtiClipFiberGroupToROIs(fg_clean(jj),roi1,roi2);
                    if length(fg_clip.fibers) > 20
                        % clean clipped fiber group
                        [~, keep]   = AFQ_removeFiberOutliers(fg_clip,afq.params.maxDist,afq.params.maxLen,afq.params.numberOfNodes,'mean',0);
                        % remove fibers from unclipped group that do no
                        % survive the clipping
                        fg_clean(jj).fibers =fg_clean(jj).fibers(keep);
                    end
                else
                    % clean un-clipped fiber group
                    fg_clean(jj) = AFQ_removeFiberOutliers(fg_clean(jj),afq.params.maxDist,afq.params.maxLen,afq.params.numberOfNodes,'mean',0,afq.params.cleanIter);
                end
            end
        end
        % Save cleaned fibers
        cleanFgName = fullfile(fibDir,[prefix(segName) '_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']);
        dtiWriteFiberGroup(fg_clean, cleanFgName);
        % Set the path to the fibers in the afq structure
        afq = AFQ_set(afq, 'clean fg path', 'subnum', ii, cleanFgName);
        % Convert fiber group back to a 1 cell structure for future
        % computations
        fg_classified = dtiFgArrayToFiberGroup(fg_clean, segName);
        
    elseif afq.params.cleanFibers == 1
        % If cleaning was already done then load the cleaned fiber group
        fprintf('\nFiber tract cleaning was already done for subject %s',sub_dirs{ii});
        fg_classified = AFQ_get(afq, 'cleaned fibers',ii);
        fg_classified = dtiFgArrayToFiberGroup(fg_classified, AFQ_get(afq,'cleanfgname',ii));  
    end
    
    %% Compute Tract Profiles
    
    if AFQ_get(afq,'compute profiles',ii)
        fprintf('\nComputing Tract Profiles for subject %s',sub_dirs{ii});
        % Determine how much to weight each fiber's contribution to the
        % measurement at the tract core. Higher values mean steaper falloff
        fWeight = AFQ_get(afq,'fiber weighting');
        % By default Tract Profiles of diffusion properties will always be
        % calculated
        [fa, md, rd, ad, cl, vol, TractProfile] = AFQ_ComputeTractProperties(...
                                                fg_classified, ...
                                                dt, ...
                                                afq.params.numberOfNodes, ...
                                                afq.params.clip2rois, ...
                                                sub_dirs{ii}, ...
                                                fWeight, ...
                                                afq);
        
        % Parameterize the shape of each fiber group with calculations of
        % curvature and torsion at each point and add it to the tract
        % profile
        [curv, tors, TractProfile] = AFQ_ParamaterizeTractShape(fg_classified, TractProfile);
        
        % Calculate the volume of each Tract Profile
        TractProfile = AFQ_TractProfileVolume(TractProfile);
        
        % Add values to the afq structure
        afq = AFQ_set(afq,'vals','subnum',ii,'fa',fa,'md',md,'rd',rd,...
            'ad',ad,'cl',cl,'curvature',curv,'torsion',tors,'volume',vol);
        
        % Add Tract Profiles to the afq structure
        afq = AFQ_set(afq,'tract profile','subnum',ii,TractProfile);
        
        % If any other images were supplied calculate a Tract Profile for that
        % parameter
        numimages = AFQ_get(afq, 'numimages');
        if numimages > 0;
            for jj = 1:numimages
                % Read the image file
                image = niftiRead(afq.files.images(jj).path{ii});
                % Check image header
                if ~all(image.qto_xyz(:) == image.sto_xyz(:))
                   image = niftiCheckQto(image);
                end  
                % Resample image to match dwi resolution if desired
                if AFQ_get(afq,'imresample')
                    image = mrAnatResampleToNifti(image, fullfile(afq.sub_dirs{ii},'bin','b0.nii.gz'),[],[7 7 7 0 0 0]);
                end
                % Compute a Tract Profile for that image
                imagevals = AFQ_ComputeTractProperties(fg_classified, image, afq.params.numberOfNodes, afq.params.clip2rois, sub_dirs{ii}, fWeight, afq);
                % Add values to the afq structure
                afq = AFQ_set(afq,'vals','subnum',ii,afq.files.images(jj).name, imagevals);
                clear imagevals
            end
        end
    else
        fprintf('\nTract Profiles already computed for subject %s',sub_dirs{ii});
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
    
    % clear the files that were computed for this subject
    clear fg fg_classified TractProfile
end

%% Compute Control Group Norms

% Check if norms should be computed
if AFQ_get(afq,'computenorms')
    % If no control group was entered then norms will only contain nans.
    [norms, patient_data, control_data, afq] = AFQ_ComputeNorms(afq);
end

%% Identify Patients With Abnormal Diffusion Measurements

property = 'FA';
% Only run AFQ_ComparePatientsToNorms if norms were computed for the
% property of interest for each tract
if AFQ_get(afq,'computenorms') == 0 || sum(isnan(eval(['norms.mean' property '(1,:)']))) == 20
    fprintf('\nnorms are empty. skipping AFQ_ComparePatientsToNorms\n')
    % If there are no norms than we can not identify which subjects are
    % abnormal.  Hence these variables will be set to nan.
    abn       = nan(length(sub_dirs),1);
    abnTracts = nan(length(sub_dirs),20);
elseif AFQ_get(afq,'number of patients') >=1
    [abn, abnTracts] = AFQ_ComparePatientsToNorms(patient_data, norms, afq.params.cutoff, property);
else
    abn = nan;
    abnTracts = nan;
end
%% Plot Abnormal Patients Against Control Population

% Only plot if norms were computed for each tract
if AFQ_get(afq,'computenorms') == 0 || sum(isnan(eval(['norms.mean' property '(1,:)']))) == 20
    fprintf('\nnorms are empty. Skipping AFQ_plot\n')
elseif ~exist('abn','var') || sum(isnan(abn))==1
    fprintf('\nNo patients. Skipping AFQ_plot\n')
elseif AFQ_get(afq,'showfigs');
    % percentiles to define normal range
    ci = afq.params.cutoff;
    % loop over tracts and plot abnormal subjects
    for jj = 1:20
        % Find subjects that show an abnormality on tract jj
        sub_nums = find(abnTracts(:,jj));
        % Generate a structure for a legend
        L = {};
        for ii = 1:length(sub_nums)
            L{ii} = num2str(sub_nums(ii));
        end
        if ~isempty(sub_nums)
            AFQ_plot(norms, patient_data,'individual','ci',ci,'subjects',sub_nums,'tracts',jj,'legend',L);
        end
        % AFQ_PlotResults(patient_data, norms, abn, afq.params.cutoff,property, afq.params.numberOfNodes, afq.params.outdir, afq.params.savefigs);
    end
end

%% Plot group means for the patients and the controls

% Only plot if there is data for patients and controls
if sum(sub_group == 1) > 2 && sum(sub_group == 0) > 2 && AFQ_get(afq,'showfigs')
    AFQ_plot('Patients', patient_data, 'Controls', control_data, 'group');
elseif sum(sub_group == 1) <= 2 && sum(sub_group == 0) <= 2
    fprintf('\nNot enough subjects for a group comparison\n')
end
