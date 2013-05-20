%% Run AFQ analysis for 2 groups (patients and controls)

% Get the path to the AFQ directories
[AFQbase AFQdata] = AFQ_directories;
% Create a cell array where each cell is the path to a data directory
sub_dirs = {[AFQdata '/patient_01/dti30'], [AFQdata '/patient_02/dti30']...
    [AFQdata '/patient_03/dti30'], [AFQdata '/control_01/dti30']...
    [AFQdata '/control_02/dti30'], [AFQdata '/control_03/dti30']};
% Create a vector of 0s and 1s defining who is a patient and a control
sub_group = [1, 1, 1, 0, 0, 0];
% Run AFQ in test mode to save time and do not show figures. No inputs are
% needed to run AFQ with the default settings. AFQ_Create builds the afq
% structure. This will also be done automatically by AFQ_run if the user
% does not wish to modify any parameters
afq = AFQ_Create('run_mode','test', 'sub_dirs', sub_dirs, 'sub_group',...
    sub_group, 'showfigs',false);
% Run AFQ to generate the fiber tracts
[afq, patient_data, control_data, norms, abn, abnTracts] = AFQ_run(sub_dirs, sub_group, afq);

%% Run a t-test to compare FA along each fiber tract for patients vs. controls

% Loop over all 20 fiber groups
for jj = 1:20
    % Run an independent samples t-test comparing FA values between the
    % groups at each point along each tract
    [h(jj,:),p(jj,:),~,Tstats(jj)] = ttest2(afq.patient_data(jj).FA,afq.control_data(jj).FA);
end

%% Make Tract Profiles of T statistics

% Load the cleaned segmented fibers for the first control subject
fg = dtiReadFibers(fullfile(sub_dirs{3},'fibers','MoriGroups_clean_D5_L4.mat'));
% Load the subject's dt6 file
dt = dtiLoadDt6(fullfile(sub_dirs{3},'dt6.mat'));
% Compute Tract Profiles with 100 nodes
numNodes = 100;
[fa, md, rd, ad, cl, volume, TractProfile] = AFQ_ComputeTractProperties(fg,dt,numNodes);
% Add the pvalues and T statistics from the group comparison to the tract
% profile. This same code could be used to add correlation coeficients or
% other statistics
for jj = 1:20
    TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','pval',p(jj,:));
    TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','Tstat',Tstats(jj).tstat);
end

%% Render The tract Profiles of T statistics

% The hot colormap is a good one because it will allow us to make regions
% of the profile where p>0.05 black.
cmap = 'hot';
% Set the color range to black out non significant regions of the tract. We
% will render the T-stat profile such that T statistics greater than 5 will
% be white and T statistics less than 1 will be black. Obviously a T value
% of 1 is not considered significant however in this example we only have 3
% subjects in each group hence we use a low value just to demonstrate the
% proceedure.
crange = [1 4];
% Set the number of fibers to render. More fibers takes longer
numfibers = 200;
% Render the left corticospinal tract (fibers colored light blue) with a
% Tract Profile of T statistics. Each fiber will have a 1mm radius and the
% tract profile will have a 6mm radius.
AFQ_RenderFibers(fg(3),'color',[.8 .8 1],'tractprofile',TractProfile(3),...
    'val','Tstat','numfibers',numfibers,'cmap',cmap,'crange',crange,...
    'radius',[1 6]);
% Add the defining ROIs to the rendering.
[roi1 roi2] = AFQ_LoadROIs(3,sub_dirs{3});
% The rois will be rendered in dark blue
AFQ_RenderRoi(roi1,[0 0 .7]);
AFQ_RenderRoi(roi2,[0 0 .7]);
% Add the slice x=-15 from the subject's b=0 image
b0 = readFileNifti(dt.files.b0);
AFQ_AddImageTo3dPlot(b0,[-15,0,0]);
