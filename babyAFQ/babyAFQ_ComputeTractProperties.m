function [fa, md, rd, ad, cl, volume, TractProfile,SuperFibersGroup,fgResampled] = babyAFQ_ComputeTractProperties(fg_classified,dt,numberOfNodes,clip2rois,subDir, roiDir, fgnum, weighting, afq)
% Compute diffusion properties along the trajectory of the fiber groups.
% The function returns vectors of diffusion properties and TractProfile
% structure for each fiber group
%
% [fa md rd ad cl volume TractProfile] = AFQ_ComputeTractProperties(fg_classified,dt,[numberOfNodes=30],[clip2rois=1],[subDir], [weighting = 1], afq)
%
% The input fiber group can either be a single fiber group or an array of
% fiber groups (See fg2Array). Each fiber group is clipped to its 2
% defining ROIs (if clip2rois == true). The fibers are sampled to a defined
% number of nodes (sampled at equidistant locations on each fiber) and
% flipped so that they start and end in an equivalent location (if
% clip2rois==false then we assume that fibers have already start and end in
% equivalent regions), then diffusion properties are calculated at each
% node as a weighted sum of the diffusion properties of each fiber at that
% node. The weights are determined based on each fibers gaussian distance
% from the core (mean) of the tract. Other quantitative maps can also be
% input into this function and Tract Profiles of these measures will be
% returned.
%
% Input arguments:
%  fg_classified  = Fiber group classified into 20 subgroups that
%                   correspond to each of the 20 tracts output by
%                   AFQ_SegmentFiberGroups
%
%  dt             = dt6 structure or an image file.
%
%  numberOfNodes  = The number of nodes to sample each fiber to.
%
%  clip2rois      = if clip2rois is set to 1 (default) then only the
%                   central portion of the fiber group spanning between
%                   the 2 defining ROIs is analyzed.
%                   Otherwise the full trajectory is analyzed.
%
%  subDir         = The directory containing the subjects dt6 file. This is
%                   needed to find the subject's ROIs if clip2rois==1. By
%                   defualt the subject's directory will be searched for
%                   based on the path given in the dt6 file (dt.dataFile).
% weighting       = By default each fiber's contribution to the
%                   measurement at the fiber tract core is weighted by its
%                   gaussian distance from the core. This can also be
%                   interpreted as the probability that it is a member of
%                   the fiber group.  weighting defines the steepness of
%                   the falloff of the weights assigned to each fiber where
%                   a value of 0 means that there is no weighting and each
%                   fiber contributes equally, a value of 1 means that we
%                   use gaussian weighting (default) and values >1 mean
%                   that the weights fall off steaper as the fibers deviate
%                   further from the core. Higher weighting values mean the
%                   core is weighted more heavily compare to the periphery
%                   of the tract.
%
% Outputs:
% fa               = vector of Fractional anisotropy along fiber core.
% md               = vector of Mean Diffusivity values along fiber core.
% rd               = vector of Radial Diffusivity values along fiber core.
% ad               = vector of Axial Diffusivity values along fiber core.
% cl               = vector of linearity values along fiber core.
% volume           = vector of volume estimates along the tract.
% TractProfile     = TractProfile structure (see AFQ_CreateTractProfile).
%                    Contains fiber group core, fiber covariance, and
%                    fiber tract properties at each node.
%
%  Example:
%
%   AFQdata = '/home/jyeatman/matlab/svn/vistadata/AFQ'
%   dt = dtiLoadDt6(fullfile(AFQdata,'subj2','dt6.mat'))
%   fg_classified = AFQ_SegmentFiberGroups(dt)
%   numberOfNodes = 30; clip2rois = 1; subDir = fullfile(AFQdata,'subj2');
%   [fa md rd ad cl SuperFibersGroup] =
%   AFQ_ComputeTractProperties(fg_classified,dt,numberOfNodes,clip2rois,subDir)
%
% Copyright Stanford Vista Team 2011
% Written by Jason D. Yeatman

% check the format of the fiber group. Convert to an array of fiber groups
% if necessary
if isstruct(fg_classified) && isfield(fg_classified,'subgroupNames')...
        && isfield(fg_classified,'subgroup')
    fg_classified = fg2Array(fg_classified);
end
% if number of nodes doesn't exist then set to 30
if ~exist('numberOfNodes','var') || isempty(numberOfNodes)
    numberOfNodes=30;
end
if ~exist('clip2rois','var') || isempty(clip2rois)
    clip2rois=1;
end
% check if the dt input is a dt6 structure or a nifti image
if isfield(dt,'qto_xyz')
    valname = 'image';
else
    valname = 'dt';
end
% if no subject directory was defined then define based on the dt6 file
if ~exist('subDir','var') || isempty(subDir) && strcmp(valname, 'dt')
    subDir = fileparts(dt.dataFile);
end
% Check if a weighting factor was defined otherwise set to gaussian
% weighting
if ~exist('weighting','var') || isempty(weighting)
    weighting = 1;
end
% Number of fiber groups
numfg = length(fg_classified);
% Create Tract Profile structure
TractProfile = AFQ_CreateTractProfile;
% Pre allocate data arrays
fa=nan(numberOfNodes,numfg);
md=nan(numberOfNodes,numfg);
rd=nan(numberOfNodes,numfg);
ad=nan(numberOfNodes,numfg);
cl=nan(numberOfNodes,numfg);
cs=nan(numberOfNodes,numfg);
cp=nan(numberOfNodes,numfg);
volume=nan(numberOfNodes,numfg);
% loop over the fiber groups
for jj=[1:6 9:26]
    % There are 20 fiber groups saved within the same structure. We will
    % loop over these fiber groups, allocate them to a tmp variable
    % calculate what we need and move on
    fgtmp = fg_classified(jj);
    
    % This is what we would do with the old fiber group structure
    % fgtmp=dtiNewFiberGroup(fg_classified.subgroupNames(jj).subgroupName);
    %  fgtmp.fibers=fg_classified.fibers(fg_classified.subgroup==jj);
    
    % Figure out the fiber group number if an afq structure was passed in
        if exist('afq','var') && ~isempty(afq)
            fgnames = AFQ_get(afq,'fgnames');
            fgnum = find(strcmpi(fgtmp.name, fgnames));
            if isempty(fgnum)
                fgnum = jj;
                fprintf('\n%s does not match any fiber group name in the afq structure.', fgtmp.name);
                fprintf('\n Assuming that it is equivalent to %s',fgnames{jj});
            end
        else
            fgnum = jj;
        end
   
    % clip the fiber group to the portion spanning between the two ROIs if
    % desired
%     if clip2rois==1 && exist('afq','var') && ~isempty(afq)
%         [roi1, roi2] = babyAFQ_LoadROIs(fgnum,subDir,afq,roiDir);
%     elseif clip2rois==1

if jj~=23 && jj~=24
        [roi1, roi2] = babyAFQ_LoadROIsReoriented(fgnum,subDir,[],roiDir);
        [fgtmp]=AFQ_ReorientFibers(fgtmp,roi1,roi2);
        
        if clip2rois~=1
        roi1=[]; roi2=[];
        end
end

%     end
    % Some subjects may be missing certain fiber groups due to lesions,
    % or noise in the data.  If a fiber group is empty we will ad NaNs to
    % the database
    if isempty(fgtmp.fibers)
        fa(:, jj)=nan;md(:, jj)=nan;rd(:, jj)=nan;ad(:, jj)=nan;
        cl(:, jj)=nan;cp(:,jj)=nan;cs(:,jj)=nan;volume(:,jj)=nan;
        TractProfile(jj) = AFQ_CreateTractProfile('name',fgtmp.name);
        continue
    end
    % Here is where we will create define the fiber group core and calculate
    % stats based on the method described in Yeatman et al 2011. Sometimes
    % there is a problem such as only having 1 fiber in the fiber group.
    % To avoid breaking the whole loop we use the try catch loop
    %try
        [fa(:, jj),md(:, jj),rd(:, jj),ad(:, jj),cl(:, jj),...
            SuperFibersGroup(jj),~,cp(:,jj),cs(:,jj),fgResampled]=...
            dtiComputeDiffusionPropertiesAlongFG(fgtmp, dt, roi1, roi2, numberOfNodes,weighting);
        % Compute tract volume at each node
        %volume(:,jj) = AFQ_TractProfileVolume(fgResampled);
        % Put mean fiber into Tract Profile structure
        TractProfile(jj) = AFQ_CreateTractProfile('name',fgtmp.name,'superfiber',SuperFibersGroup(jj));
        % Add the volume measurement
        TractProfile(jj).fiberVolume = volume(:,jj)';
        % Add planarity and sphericity measures to the tract profile. We
        % could return them as outputs from the main function but the list
        % of outputs keeps growing...
        TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','planarity',cp(:,jj)');
        TractProfile(jj) = AFQ_TractProfileSet(TractProfile(jj),'vals','sphericity',cs(:,jj)');
%     catch
%          fa(:, jj)=nan;md(:, jj)=nan;rd(:, jj)=nan;ad(:, jj)=nan;
%          cl(:, jj)=nan;cp(:,jj)=nan;cs(:,jj)=nan;volume(:,jj) = nan;
%          TractProfile(jj) = AFQ_CreateTractProfile('name',fgtmp.name);
%     end 
end


return