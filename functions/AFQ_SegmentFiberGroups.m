function [fg_classified,fg_unclassified,classification,fg] = ...
    AFQ_SegmentFiberGroups(dt6File, fg, Atlas, ...
    useRoiBasedApproach, useInterhemisphericSplit, antsInvWarp, template)
% Categorizes each fiber in a group into one of the 20 tracts defined in
% the Mori white matter atlas. 
%
%  [fg_classified, fg_unclassified, classification, fg] = ...
%      AFQ_SegmentFiberGroups(dt6File, fg, [Atlas='MNI_JHU_tracts_prob.nii.gz'], ...
%      [useRoiBasedApproach=true], [useInterhemisphericSplit=true], [antsInvWarp] , [template]);
%
%  Fibers are segmented in two steps. Fibers become candidates for a fiber
%  group if the pass through the 2 waypoint ROIs that define the
%  tracjectory of the tract. Then each fiber is compared to a fiber
%  proability map and high probability fibers are retained in the group.
%  The segmentation alogrithm is based on:
%
%  Hua K, Zhang J, Wakana S, Jiang H, Li X, Reich DS, Calabresi PA, Pekar
%  JJ, van Zijl PC, Mori S. 2008. Tract probability maps in stereotaxic
%  spaces: analyses of white matter anatomy and tract-specific
%  quantification. Neuroimage 39(1):336-47.
%
%  Zhang, W., Olivi, A., Hertig, S. J., van Zijl, P., & Mori, S. (2008).
%  Automated fiber tracking of human brain white matter using diffusion
%  tensor imaging. Neuroimage, 42(2), 771-777.
%
% Input parameters:
% dt6File                  - Either the dt6 structure or a path to the
%                            dt6.mat file
% fg                       - A file with previously tracked elsewhere fibers
%                          to be categorized.
% Atlas                    - probabilistic atlas defining probabilities for
%                          each voxel to be passed by a fiber within each of
%                          atlas fiber groups. We usually use Mori atlas
%                          supplied with fsl: MNI_JHU_tracts_prob.nii.gz.
%                          This atlas is not symmetric. For a "symmetrified"
%                          atlas use 'MNI_JHU_tracts_prob_Symmetric.nii.gz'
%                          but we strongly recommend using the original
%                          atlas.
% useInterhemisphericSplit - cut fibers crossing between hemispheres with a
%                          midsaggital plane below z=-10. This is to get
%                          rid of CST fibers that appear to cross at the
%                          level of the brainstem
% useRoiBasedApproach      - use the approach describing in Zhang (2008) Neuroimage 42.
%                           For each of the 20 Mori Groups 2 critical ROIs
%                           are computed by spatially transforming ROIs
%                           provided in templates/MNI_JHU_tracts_ROIs. A
%                           fiber becomes a candidate to being labeled as a
%                           part of a given Mori group if this fiber
%                           "passes through" both critical ROIs for that
%                           Mori group. Our modification of Zhang (2008)
%                           approach: In case a single fiber is a candidate
%                           for >1 Mori group, respective cumulative
%                           probabilities are computed with probabilistic
%                           Mori atlas, then compared.
%                           useRoiBasedApproach can take the following
%                           values: (1) 'false' (to not use the approach);
%                           (2) 'true'(use the approach; the minimal
%                           distance from a fiber to ROI to count as "a
%                           fiber  is crossing the ROI" minDist=2mm); (3) a
%                           scalar value for minDist in mm; (4) a vector
%                           with the first value being minDist and the
%                           second value being a flag 1/0 for whether ROIs
%                           should be recomputed (and overwritten) for this
%                           subject (1, by default). This is useful because
%                           sometimes if you  are rerunning
%                           dtiFindMoriTracts with different parameters,
%                           you do not need to recompute ROIs.  E.g., to
%                           avoid recomputing  ROIs and use minDist of 4mm
%                           one would pass [useRoiBasedApproach=[4 0]];
%  antsInvWarp              - Spatial normalization computed with ANTS. If
%                           a path to a precomputed ANTS warp is passed in
%                           then it will be used to transform the MNI ROIs
%                           to native space
% template                  - Path to an MNI template
%
% Output parameters:
% fg_ classified  - fibers structure containing all fibers assigned to
%                   one of Mori Groups. Respective group labeles are stored
%                   in fg.subgroups field.
% fg_unclassified - fiber structure containing the rest of the (not Mori) fibers.
% classification  - This variable gives the fiber group names and the group
%                   fiber group number for each fiber in the input group
%                   fg.  classification is a structure with two fields.
%                   classification.names is a cell array where each cell is
%                   the name of that fiber group. For example
%                   classification.names{3} = 'Corticospinal tract L'.
%                   classification.index is a vector that defines which
%                   group number each fiber in the origional fiber group
%                   was assigned to. For example
%                   classification.index(150)=3 means that fg.fibers(150)
%                   is part of the corticospinal tract fiber group.  The
%                   values in classification may not match the origional
%                   fiber group because of pre-processing.  However they
%                   will match the output fg which is the origional group
%                   with preprocessing.
% fg              - This is the origional pre-segmented fiber group.  It
%                   may differ slightly from the input due to preprocessing
%                   (eg splitting fibers that cross at the pons, removing 
%                   fibers that are too short)
%                    
% Example:
%    AFQdata = '/home/jyeatman/matlab/svn/vistadata/AFQ';
%    dt6File = fullfile(AFQdata, 'subj2', 'dt6.mat');
%    fg      = AFQ_SegmentFiberGroups(dt6File);
%
% See also: dtiSplitInterhemisphericFibers
%
% (c) Vistalab

%% Check arguments

if ~exist('useInterhemisphericSplit', 'var') || isempty(useInterhemisphericSplit)
    useInterhemisphericSplit=true;
end

if ~exist('useRoiBasedApproach', 'var')|| isempty(useRoiBasedApproach)
    useRoiBasedApproach=true;
end
% Decide on the recomputeROI condition
if useRoiBasedApproach==false,         recomputeROIs = false;
elseif length(useRoiBasedApproach)<2,  recomputeROIs = 1;
else                                   recomputeROIs = useRoiBasedApproach(2);
end
if recomputeROIs
    display('You chose to recompute ROIs');
end
%E.g., to avoid recomputing  ROIs and use minDist of 4mm one would pass [useRoiBasedApproach=[4 0]];
if isnumeric(useRoiBasedApproach)
    minDist=useRoiBasedApproach(1);
    useRoiBasedApproach='true';
else
    minDist=2; %defualt is .89;
end
display(['Fibers that get as close to the ROIs as ' num2str(minDist) 'mm will become candidates for the Mori Groups']);
% This is left as an option in case an updated atlas is released
if(~exist('Atlas','var') || isempty(Atlas))
    %Default scenario: use original Mori Atlas
    Atlas='MNI_JHU_tracts_prob.nii.gz';
end

if ~exist('template', 'var') || isempty(template)
    tdir = fullfile(AFQ_directories,'templates','mni_icbm152_nlin_asym_09a_nifti');
    template = fullfile(tdir,'mni_icbm152_t2_tal_nlin_asym_09a.nii');
end
%% Read the data
% Load the dt6 file
if ischar(dt6File)
    dt = dtiLoadDt6(dt6File);
    baseDir = fileparts(dt6File); 
else
    dt = dt6File;
    baseDir = fileparts(dt.dataFile);
    dt6File = dt.dataFile;
end
% Track wholebrain fiber group if one was not passed in
if ~exist('fg','var') || isempty(fg)
    fg = AFQ_WholebrainTractography(dt);
end
% Load fiber group - Can be filename or the data
if ischar(fg), fg = dtiLoadFiberGroup(fg); end
% Initialize spm defualts for normalization
spm_defaults; global defaults; params = defaults.normalise.estimate;
% spm_get_defaults - For SPM8

%% Spatially normalize diffusion data with the MNI (ICBM) template
% Rescale image valueds to get better gary/white/CSF contrast
alignIm = mrAnatHistogramClip(double(dt.b0),0.3,0.99);
% Compute normalization
[sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(alignIm, dt.xformToAcpc, template, params);
% check the normalization
mm = diag(chol(Vtemplate.mat(1:3,1:3)'*Vtemplate.mat(1:3,1:3)))';
bb = mrAnatXformCoords(Vtemplate.mat,[1 1 1; Vtemplate.dim]);
alignIm_sn = mrAnatResliceSpm(alignIm, sn, bb, [2 2 2], [1 1 1 0 0 0], 0);
tIm = mrAnatResliceSpm(double(Vtemplate.dat), inv(Vtemplate.mat), bb, [2 2 2], [1 1 1 0 0 0], 0);
im(:,:,:,1) = uint8(tIm);
im(:,:,:,2) = uint8(round(clip(alignIm_sn)*255));
im(:,:,:,3) = im(:,:,:,2);
% Save an image of the spatially normalized data showing the quality of the
% alignment
imwrite(makeMontage(im),fullfile(baseDir, 'SpatialNormalization.png'));

%% Compute fiber group probabilities using the atlas proceedure of Hua.2008
% Load the Mori atlas maps these are saved in nifti images
tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');
moriTracts = readFileNifti(fullfile(tdir, Atlas));
% 15 is a subregion of 19 and 16 a subregion of 20. To better separate them,
% we subtract 19 from 15 and 20 from 16.
moriTracts.data(:,:,:,15) = moriTracts.data(:,:,:,15)-moriTracts.data(:,:,:,19);
moriTracts.data(:,:,:,16) = moriTracts.data(:,:,:,16)-moriTracts.data(:,:,:,20);
% Load the fiber group labels
%labels = readTab(fullfile(tdir,'MNI_JHU_tracts_prob.txt'),',',false);
%labels = labels(1:20,2);
labels = {'Left Thalamic Radiation','Right Thalamic Radiation','Left Corticospinal','Right Corticospinal', 'Left Cingulum Cingulate', 'Right Cingulum Cingulate'...
    'Left Cingulum Hippocampus','Right Cingulum Hippocampus', 'Callosum Forceps Major', 'Callosum Forceps Minor'...
    'Left IFOF','Right IFOF','Left ILF','Right ILF','Left SLF','Right SLF','Left Uncinate','Right Uncinate','Left Arcuate','Right Arcuate'};

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% If you wanted to inverse-normalize the maps to this subject's brain:
% invDef.outMat = moriTracts.qto_ijk;
% bb = mrAnatXformCoords(dt.xformToAcpc,[1 1 1; size(dt.b0)]);
% tprob = mrAnatResliceSpm(tprob, invDef, bb, dt.mmPerVoxel, [1 1 1 0 0 0]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Cut the fibers below the acpc plane, to disentangle CST & ATL crossing 
% at the level of the pons
if useInterhemisphericSplit
    fgname  = fg.name;
    fg      = dtiSplitInterhemisphericFibers(fg, dt, -10);
    fg.name = fgname;
end
% Throw out fibers that are too short to span between two ROIs
if sum(cellfun(@length, fg.fibers)<5)~=0
    if isfield(fg, 'subgroup')&&~isempty(fg.subgroup)
        fg.subgroup(cellfun(@length, fg.fibers)<5)=[];
    end
    if isfield(fg, 'seeds')&&~isempty(fg.seeds)
        fg.seeds(cellfun(@length, fg.fibers)<5, :)=[];
    end
    fg.fibers(cellfun(@length, fg.fibers)<5)=[];
end

% Warp the fibers to the MNI standard space so they can be compared to the
% template
fg_sn = dtiXformFiberCoords(fg, invDef);

% moriTracts.data is a an XxYxZx20 array contianing the 20 Mori probability
% atlases (range is 0-100 where 100 represents p(1)).
sz = size(moriTracts.data);
% fg_sn fiber coords are in MNI space- now convert them to atlas space by
% applying the affine xform from the atlas NIFTI header. Since the atlas is
% already in MNI space, this transform will just account for any
% translation and scale differences between the atlas maps and the MNI
% template used to compute our sn.
fgCoords = mrAnatXformCoords(moriTracts.qto_ijk, horzcat(fg_sn.fibers{:}));
clear fg_sn;   % what we need from fg_sn is now stored in fgCoords
fgLen = cellfun('size',fg.fibers,2);

% Now loop over the 20 atlases and get the atlas probability score for each
% fiber point. We collapse the scores across all points in a fiber by taking
% the mean. Below, we will use these 20 mean scores to categorize the fibers.
% TO DO: consider doing something more sophisticated than taking the mean.
fp = zeros(sz(4),numel(fg.fibers));
for(ii=1:sz(4))
    % Get the Mori atlas score for each point in the fibers using
    % trilinear interpolation.
    p = myCinterp3(double(moriTracts.data(:,:,:,ii))/100, sz([1,2]), sz(3), fgCoords(:,[2,1,3]));
    % The previous line interpolated one giant array with all fiber points
    % concatenated. The next loop will separate the coordinates back into
    % fibers and take the mean score for the points within each fiber.
    fiberCoord = 1;
    for(jj=1:numel(fg.fibers))
        fp(ii,jj) = nanmean(p([fiberCoord:fiberCoord+fgLen(jj)-1]));
        fiberCoord = fiberCoord+fgLen(jj);
    end
end
clear p fgCoords;

%% Find fibers that pass through both waypoint ROIs eg. Wakana 2007
%Warp Mori ROIs to individual space; collect candidates for each fiber
%group based on protocol of 2 or > ROIs a fiber should travel through. The
%following ROIs are saved within
%trunk/mrDiffusion/templates/MNI_JHU_tracts_ROIs folder and are created
%using MNI template as described in Wakana et al.(2007) Neuroimage 36 with
%a single modification: For SLFt Roi2, they recommend drawing the ROI at
%the AC level, whereas we use a slice just inferior of CC splenium. The
%reason for this modification is that Wakana et al. ACPC aligned images
%appear different from MNI images (the latter we use for defininng ROIs).
%If defining SLFt-Roi2 on a slice actually at the AC level (althought
%highly consistently across human raters), many SLFt fibers were not
%correctly labeled as they extend laterally into temporal lobe just above
%the aforementioned ROI plane.
if useRoiBasedApproach 
    % A 20x2 cell array containing the names of both waypoint ROIs for each
    % of 20 fiber groups
    moriRois={'ATR_roi1_L.nii.gz',  'ATR_roi2_L.nii.gz'; 'ATR_roi1_R.nii.gz', 'ATR_roi2_R.nii.gz'; ...
        'CST_roi1_L.nii.gz', 'CST_roi2_L.nii.gz'; 'CST_roi1_R.nii.gz',  'CST_roi2_R.nii.gz'; ...
        'CGC_roi1_L.nii.gz', 'CGC_roi2_L.nii.gz'; 'CGC_roi1_R.nii.gz', 'CGC_roi2_R.nii.gz'; ...
        'HCC_roi1_L.nii.gz', 'HCC_roi2_L.nii.gz'; 'HCC_roi1_R.nii.gz', 'HCC_roi2_R.nii.gz';...
        'FP_R.nii.gz', 'FP_L.nii.gz'; ...
        'FA_L.nii.gz', 'FA_R.nii.gz'; ...
        'IFO_roi1_L.nii.gz', 'IFO_roi2_L.nii.gz'; 'IFO_roi1_R.nii.gz', 'IFO_roi2_R.nii.gz'; ...
        'ILF_roi1_L.nii.gz', 'ILF_roi2_L.nii.gz'; 'ILF_roi1_R.nii.gz', 'ILF_roi2_R.nii.gz'; ...
        'SLF_roi1_L.nii.gz', 'SLF_roi2_L.nii.gz'; 'SLF_roi1_R.nii.gz', 'SLF_roi2_R.nii.gz'; ...
        'UNC_roi1_L.nii.gz', 'UNC_roi2_L.nii.gz'; 'UNC_roi1_R.nii.gz', 'UNC_roi2_R.nii.gz'; ...
        'SLF_roi1_L.nii.gz', 'SLFt_roi2_L.nii.gz'; 'SLF_roi1_R.nii.gz', 'SLFt_roi2_R.nii.gz'};
    % Make an ROI for the mid saggital plane
    midSaggitalRoi= dtiRoiMakePlane([0, dt.bb(1, 2), dt.bb(1, 3); 0 , dt.bb(2, 2) , dt.bb(2, 3)], 'midsaggital', 'g');
    keep1=zeros(length(fg.fibers), size(moriRois, 1)); keep2=zeros(length(fg.fibers), size(moriRois, 1));
    % Find fibers that cross mid saggital plane
    [fgOut, contentiousFibers, InterHemisphericFibers] = dtiIntersectFibersWithRoi([], 'not', [], midSaggitalRoi, fg); 
    %NOTICE: ~keep3 (not "keep3") will mark fibers that DO NOT cross
    %midSaggitalRoi.
    keep3=repmat(InterHemisphericFibers, [1 size(moriRois, 1)]);
    fgCopy=fg; fgCopy.subgroup=[];
    for roiID=1:size(moriRois, 1)
        % Load the nifit image containing ROI-1 in MNI space
        ROI_img_file=fullfile(tdir, 'MNI_JHU_tracts_ROIs',  [moriRois{roiID, 1}]);
        % Transform ROI-1 to an individuals native space
        if recomputeROIs
            % Default is to use the spm normalization unless a superior
            % ANTS normalization was passed in
            if exist('antsInvWarp','var') && ~isempty(antsInvWarp)
                outfile = fullfile(fileparts(dt6File),'ROIs',moriRois{roiID, 1});
                roi1(roiID) = ANTS_CreateRoiFromMniNifti(ROI_img_file, antsInvWarp, [], outfile);
            else
                [RoiFileName, invDef, roi1(roiID)]=dtiCreateRoiFromMniNifti(dt6File, ROI_img_file, invDef, true);
            end
        else
            RoiFileName=fullfile(fileparts(dt6File), 'ROIs',  [prefix(prefix(ROI_img_file, 'short'), 'short') '.mat']);
            roi1(roiID) = dtiReadRoi(RoiFileName);  
        end
        % Find fibers that intersect the ROI
        [fgOut,contentiousFibers, keep1(:, roiID)] = dtiIntersectFibersWithRoi([], 'and', minDist, roi1(roiID), fg);
        keepID1=find(keep1(:, roiID));
        % Load the nifit image containing ROI-2 in MNI space
        ROI_img_file=fullfile(tdir, 'MNI_JHU_tracts_ROIs',  [moriRois{roiID, 2}]);
        % Transform ROI-2 to an individuals native space
        if recomputeROIs
            [RoiFileName, invDef, roi]=dtiCreateRoiFromMniNifti(dt6File, ROI_img_file, invDef, true);
            % Default is to use the spm normalization unless a superior
            % ANTS normalization was passed in
            if exist('antsInvWarp','var') && ~isempty(antsInvWarp)
                outfile = fullfile(fileparts(dt6File),'ROIs',moriRois{roiID, 2});
                roi2(roiID) = ANTS_CreateRoiFromMniNifti(ROI_img_file, antsInvWarp, [], outfile);
            else
                [RoiFileName, invDef, roi2(roiID)]=dtiCreateRoiFromMniNifti(dt6File, ROI_img_file, invDef, true);
            end   
        else
            RoiFileName=fullfile(fileparts(dt6File), 'ROIs',  [prefix(prefix(ROI_img_file, 'short'), 'short') '.mat']);
            roi2(roiID) = dtiReadRoi(RoiFileName);
        end
        %To speed up the function, we intersect with the second ROI not all the
        %fibers, but only those that passed first ROI.
        fgCopy.fibers=fg.fibers(keepID1(keepID1>0));
        [a,b, keep2given1] = dtiIntersectFibersWithRoi([], 'and', minDist, roi2(roiID), fgCopy);
        keep2(keepID1(keep2given1), roiID)=true; 
    end
    clear fgOut contentiousFibers keepID
    %Note: forceps major and minor should NOT have interhemipsheric fibers
    % excluded
    keep3(:, 9:10)=keep3(:, 9:10).*0;
    % fp is the variable containing each fibers score for matching the
    % atlas. We will set each fibers score to 0 if it does not pass through
    % the necessary ROIs
    fp(~(keep1'&keep2'&~keep3'))=0;
    %Also note: Tracts that cross through slf_t rois should be automatically
    %classified as slf_t, without considering their probs.
    fp(19, (keep1(:, 19)'&keep2(:, 19)'&~keep3(:, 19)'))=max(fp(:));
    fp(20, (keep1(:, 20)'&keep2(:, 20)'&~keep3(:, 20)'))=max(fp(:));
end

%% Eliminate fibers that don't match any of the atlases very well
% We have a set of atlas scores for each each fiber. To categorize the
% fibers, we will find the atlas with the highest score (using 'sort').
[atlasScore,atlasInd] = sort(fp,1,'descend');
% Eliminate fibers that don't match any of the atlases very well:
unclassified=atlasScore(1,:)==0;
goodEnough = atlasScore(1,:)~=0;
curAtlasFibers = cell(1,sz(4));
for ii=1:sz(4)
    curAtlasFibers{ii} = find(atlasInd(1,:)==ii & goodEnough);
end

% We now have a cell array (curAtlasFibers) that contains 20 arrays, each
% listing the fiber indices for the corresponding atlas group. E.g.,
% curAtlasFibers{3} is a list of indices into fg.fibers that specify the
% fibers belonging to group 3.

%Create a FG for unclassified fibers
fg_unclassified=fg;
fg_unclassified.name=[fg.name ' not Mori Groups'];
fg_unclassified.fibers = fg.fibers(unclassified); %prepare fg output
if ~isempty(fg.seeds)
    fg_unclassified.seeds = fg.seeds(unclassified);
end
fg_unclassified.subgroup=zeros(size(fg.fibers(unclassified)))'+(1+sz(4));
fg_unclassified.subgroupNames(1)=struct('subgroupIndex', 1+sz(4), 'subgroupName', 'NotMori');

% Create a fiber group for classified fibers
fg_classified = fg;
% Modify fg.fibers to discard the fibers that didn't make it into any of
% the atlas groups:
fg_classified.name=[fg.name ' Mori Groups'];
fg_classified.fibers = fg.fibers([curAtlasFibers{:}]);
if ~isempty(fg.seeds)
    fg_classified.seeds = fg.seeds([curAtlasFibers{:}],:);
end

% We changed the size of fg_classified.fibers by discarding the
% uncategorized fibers, so we need to create a new array to categorize the
% fibers. This time we make an array with one entry corresponding to each
% fiber, with integer values indicating to which atlas group the
% corresponding fiber belongs.
fg_classified.subgroup = zeros(1,numel(fg_classified.fibers));
curInd = 1;
for(ii=1:numel(curAtlasFibers))
    fg_classified.subgroup(curInd:curInd+numel(curAtlasFibers{ii})-1) = ii;
    %fghandle.subgroup(curInd:curInd+numel(curAtlasFibers{ii})-1) = ii;
    curInd = curInd+numel(curAtlasFibers{ii});
    %Save labels for the fiber subgroups within the file
    fg_classified.subgroupNames(ii)=struct('subgroupIndex', ii, 'subgroupName', labels(ii));
    %fghandle.subgroupNames(ii)=struct('subgroupIndex', ii, 'subgroupName', labels(ii));
end
% The number of fibers in the origional preprocessed fiber group (fg) is
% equalt to the number of fibers in fg_classified + fg_unclassified. To
% create a new structure denoting the classification of each fiber group in
% the origional structure first preallocate an array the length of the
% origional (preprocessed) fiber group.
fiberIndex = zeros(length(fg.fibers),1);

% Populate the structure denoting the fiber group number that each fiber in
% the origional wholebrain group was assigned to
for ii = 1 : length(curAtlasFibers)
    fiberIndex(curAtlasFibers{ii}) = ii;
    names{ii} = fg_classified.subgroupNames(ii).subgroupName;
end

% Create a structure with fiber indices and group names.
classification.index = fiberIndex;
classification.names = names;

% Convert the fiber group into an array of fiber groups
fg_classified = fg2Array(fg_classified);

% Flip fibers so that each fiber in a fiber group passes through roi1
% before roi2
for ii = 1:size(moriRois, 1)
    fg_classified(ii) = AFQ_ReorientFibers(fg_classified(ii),roi1(ii),roi2(ii));
end

return;






