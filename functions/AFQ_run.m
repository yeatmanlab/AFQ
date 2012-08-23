function [patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq)
% Run AFQ analysis on a set of subjects.
%
% [patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, [afq])
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
%   % with the default settings 
%   afq = AFQ_Create('run_mode','test'); 
%   [patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq)
%
% Copyright Stanford Vista Team 2011 Written by Jason D. Yeatman, Brian A.
% Wandell and Robert F. Dougherty
%

%% Check Inputs
if notDefined('sub_dirs'), error('No subject directories');  end
if ischar(sub_dirs,1), sub_dirs = cellstr(sub_dirs); end
if ~exist('sub_group', 'var') || isempty(sub_group), error('Must define subject group'); end
if length(sub_group) ~= size(sub_dirs,1)
    error('Mis-match between subject group description and subject data directories');
end
if ~exist('afq','var') || isempty(afq)
    afq = AFQ_Create; %if no parameters are defined use the defualts
end
if isempty(afq.sub_group)
    afq = AFQ_set(afq,'sub_group',sub_group);
end
%%  Loop over every subject
for ii=1:length(sub_dirs)
    % Define the current subject to process
    afq = AFQ_set(afq,'current subject',ii);
    %% Preprocess Data
    
    if ~exist(fullfile(sub_dirs{ii},'dt6.mat'),'file')
        error('AFQ preprocessing has not yet been implemented')
    end
    
    %% Perform Whole Brain Streamlines Tractography
    % Load Dt6 File
    dtFile = fullfile(sub_dirs{ii},'dt6.mat');
    dt     = dtiLoadDt6(dtFile);
    %Check if whole=brain tractography has already been done
    fibDir = fullfile(sub_dirs{ii},'fibers');
    if ~exist(fibDir,'dir')
        mkdir(fibDir);
    end   
    % If some of the steps have already been done, there will be a
    % WholeBrainFG.mat file.  So use that one.
    if exist(fullfile(fibDir,'MoriGroups.mat'),'file')
        fprintf('\nFound segmented fiber groups for subject %s',sub_dirs{ii});
    elseif exist(fullfile(fibDir,'WholeBrainFG.mat'),'file')
        fg = dtiLoadFiberGroup(fullfile(fibDir,'WholeBrainFG.mat'));
    else
        % Perform whole-brain tractography
        fprintf('\nPerforming whole-brain tractograpy for subject %s\n',sub_dirs{ii});
        fg = AFQ_WholebrainTractography(dt, afq.params.run_mode);
        % Save fiber group to fibers directory
        dtiWriteFiberGroup(fg,fullfile(fibDir,'WholeBrainFG.mat'));
    end
    
    %% Segment 20 Fiber Groups
    
    if afq.params.cleanFibers == 1 && exist(fullfile(fibDir,['MoriGroups_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']),'file')
        fprintf('\nFound cleaned groups for subject %s',sub_dirs{ii});
    elseif exist(fullfile(fibDir,'MoriGroups.mat'),'file')
        fg_classified = dtiLoadFiberGroup(fullfile(fibDir,'MoriGroups.mat'));
    else
        % Segment fiber file
        fg_classified = AFQ_SegmentFiberGroups(dtFile, fg);
        % Save segmented fiber group
        dtiWriteFiberGroup(fg_classified, fullfile(fibDir,'MoriGroups.mat'));
        clear fg
    end
    clear dtFile fgFile
    
    %% Remove outliers from each fiber tract so it is a compact bundle
    
    if afq.params.cleanFibers == 1 && exist(fullfile(fibDir,['MoriGroups_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']),'file')
        fg_classified = dtiLoadFiberGroup(fullfile(fibDir,['MoriGroups_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']));
        fg_classified = dtiFgArrayToFiberGroup(fg_classified, 'MoriGroups');
    elseif afq.params.cleanFibers == 1
        % Remove all fibers that are too long and too far from the core of
        % the group.  This algorithm will constrain the fiber group to
        % something that can be reasonable represented as a 3d gaussian
        fg_clean = fg2Array(fg_classified); % fiber groups into an array   
        % remove fiber outliers
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
                    fg_clean(jj) = AFQ_removeFiberOutliers(fg_clean(jj),afq.params.maxDist,afq.params.maxLen,afq.params.numberOfNodes,'mean',0);
                end
            end
        end
        % Save cleaned fibers
        dtiWriteFiberGroup(fg_clean, fullfile(fibDir,['MoriGroups_clean_D' num2str(afq.params.maxDist) '_L'  num2str(afq.params.maxLen) '.mat']));
        fg_classified = dtiFgArrayToFiberGroup(fg_clean, 'MoriGroups');
    end
    
    %% Compute Tract Profiles
    
    % By default Tract Profiles of diffusion properties will always be
    % calculated
    [fa md rd ad cl] = AFQ_ComputeTractProperties(fg_classified, dt, afq.params.numberOfNodes, afq.params.clip2rois, sub_dirs{ii});
    
    % Take the stats that were calculated in the previous function and add
    % them to a sructure for the full sample of subjects.  Each fiber group
    % has its own cell. Within each cell there is a row for each subject with
    % numberofnodes columns
    % TODO: make object oriented and get rid of this
    for jj = 1:20
        groupFA{jj}(ii,:) = fa(:, jj);
        groupMD{jj}(ii,:) = md(:, jj);
        groupRD{jj}(ii,:) = rd(:, jj);
        groupAD{jj}(ii,:) = ad(:, jj);
        groupCL{jj}(ii,:) = cl(:, jj);
    end
    % Add values to the afq structure
    afq = AFQ_set(afq,'vals','subnum',ii,'fa',fa,'md',md,'rd',rd,'ad',ad,'cl',cl);
    
    % If any other images were supplied calculate a Tract Profile for that
    % parameter
    numimages = AFQ_get(afq, 'numimages');
    if numimages > 0;
        for jj = 1:numimages
            image = readFileNifti(afq.files.images(jj).path{ii});
            imagevals = AFQ_ComputeTractProperties(fg_classified, image, afq.params.numberOfNodes, afq.params.clip2rois, sub_dirs{ii});
            afq = AFQ_set(afq,'vals', afq.files.images(jj).name, imagevals);
            clear imagevals
        end
    end
    
end

%% Generate Control Group Norms

% If no control group was entered then norms will only contain nans.
[norms patient_data control_data afq] = AFQ_ComputeNorms(afq);

%% Identify Patients With Abnormal Diffusion Measurements

property = 'FA';
% Only run AFQ_ComparePatientsToNorms if norms were computed for the
% property of interest for each tract
if sum(isnan(eval(['norms.mean' property '(1,:)']))) == 20
    fprintf('\nnorms are empty. skipping AFQ_ComparePatientsToNorms\n')
    % If there are no norms than we can not identify which subjects are
    % abnormal.  Hence these variables will be set to nan.
    abn       = nan(length(sub_dirs),1);
    abnTracts = nan(length(sub_dirs),20);
elseif AFQ_get(afq,'number of patients') >=1
    [abn abnTracts] = AFQ_ComparePatientsToNorms(patient_data, norms, afq.params.cutoff, property);
else
    abn = nan;
    abnTracts = nan;
end
%% Plot Abnormal Patients Against Control Population

% Only plot if norms were computed for each tract
if sum(isnan(eval(['norms.mean' property '(1,:)']))) == 20
    fprintf('\nnorms are empty. Skipping AFQ_plot\n')
elseif ~exist('abn','var') || isnan(abn)
    fprintf('\nNo patients. Skipping AFQ_plot\n')
else
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
        AFQ_plot(norms, patient_data,'individual','ci',ci,'subjects',sub_nums,'tracts',jj,'legend',L)
        % AFQ_PlotResults(patient_data, norms, abn, afq.params.cutoff,property, afq.params.numberOfNodes, afq.params.outdir, afq.params.savefigs);
    end
end

%% Plot group means for the patients and the controls

% Only plot if there is data for patients and controls
if sum(sub_group == 1) > 2 && sum(sub_group == 0) > 2
    AFQ_plot('Patients', patient_data, 'Controls', control_data, 'group');
else
    fprintf('\nNot enough subjects for a group comparison\n')
end
