function [fg, keep, flipped] = AFQ_DefineFgEndpoints(fg,startpoint,endpoint,dt6Path,dCrit,invDef)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Define the cortical endpoints of a fiber group, remove fibers that do not
% terminate in the correct location and reorient each fiber in the group so
% that it has consistent start and endpoints. Cortical regions are rough
% definitions based on an MNI template
%
% [fg keep flipped] = AFQ_DefineFgEndpoints(fg,startpoint,endpoint,dt6Path,[dCrit],[invDef])
%
% Inputs:
%
% Outputs:
%
%
% Copyright Jason D. Yeatman October 2012

if ~exist('fg','var') || isempty(fg)
    error('please supply a fiber group')
elseif isfield(fg,'subgroup') && isfield(fg,'subgroupNames')
    array = 0;
    fgName = fg.name;
    % Convert to an array if necessary
    fg = fg2Array(fg);
end 

% Number of fiber groups
nFG = length(fg);

% If there are 20 fiber groups in the structure and no start or endpoints
% were defined then assume it is the moriGroups and handle apropriately
if nFG == 20 && (~exist('startpoint','var') || isempty(startpoint)) && ...
        (~exist('endpoint','var') || isempty(endpoint))
    startpoint = {'Thalamus_L' 'Thalamus_R' 'bcerebellum' 'cerebellum' ...
        'Cingulum_Post_L' 'Cingulum_Post_R' 'Hippocampus_L' 'Hippocampus_R'...
        'leftoccipital' 'leftfrontal' 'leftoccipital' 'rightoccipital' ...
        'leftoccipital' 'rightoccipital' 'leftinfparietal' 'rightinfparietal'...
        'leftanttemporal' 'rightanttemporal' 'leftfrontal' 'rightfrontal'};       
end

if ~exist('dt6Path','var') || isempty(dt6Path)
    error('please input a dt6.mat file')
elseif ischar(dt6Path)
    dt = dtiLoadDt6(dt6Path);
elseif isstruct(dt6Path)
    dt = dt6Path;
    dt6Path = dt.dataFile;
end
% Distance criteria for endpoints to match ROIs
if ~exist('dCrit','var') || isempty(dCrit)
    if isempty(endpoint)
        % If no endoint rois were defined then assume that fibers should
        % just be flipped but not removed
        dCritSq = 1000000;
    else
        dCritSq = 16;
    end
else
    dCritSq = dCrit^2;
end
% Compute normalization if it was not passed in
if ~exist('invDef','var') || isempty(invDef)
    % MNI template
    template = fullfile(fileparts(which('mrDiffusion.m')), 'templates','MNI_EPI.nii.gz');
    % Compute spatial normalization
    [~, ~, invDef] = mrAnatComputeSpmSpatialNorm(dt.b0, dt.xformToAcpc, template);
end
% Get the AFQ base directory
AFQbase = AFQ_directories;
% Template directory
tdir = fullfile(AFQbase,'templates','labelMaps');
% AAL template
AAL = fullfile(tdir,'MNI_AAL.nii.gz');

%% Reorient the fibers within each suplied fiber group
% Pre allocate variables
keep = cell(length(fg),1);
flipped = cell(length(fg),1);
for jj = 1:nFG
    % Check if the fiber group is empty
    if isempty(fg(jj).fibers)
        continue
    end
    % Get the numbers of numbers of the regions that define the desired ROI
    % in the AAL template
    if iscell(startpoint)
        Lnum1 = getLabelNumber(startpoint{jj});
    else
        Lnum1 = getLabelNumber(startpoint);
    end
    if iscell(endpoint)
        Lnum2 = getLabelNumber(endpoint{jj});
    else
        Lnum2 = getLabelNumber(endpoint);
    end
    % Make ROIs in native space from the template ROIs Load the desired
    % startpoint ROI from the AAL template and transform it
    % into the individual's native space
    [~, invDef, roiStart] = dtiCreateRoiFromMniNifti(dt6Path,AAL,invDef,0,Lnum1);
    % Only create an endpoint ROI if endpoints were defined
    if ~isempty(Lnum2)
        [~,~, roiEnd]  = dtiCreateRoiFromMniNifti(dt6Path,AAL,invDef,0,Lnum2);
    end
    % Make a function to compute the distance between fiber endpoints and
    % the desired start and endpoint
    distfun1 = @(x) nearpoints([x(:,1) x(:,end)],roiStart.coords');
    % Compute distance for each fiber to the start and end ROIs
    [~, dist1] = cellfun(distfun1,fg(jj).fibers,'UniformOutput',false);
    % Do the same for the endpoints if they were defined
    if ~isempty(Lnum2)
        distfun2 = @(x) nearpoints([x(:,1) x(:,end)],roiEnd.coords');
        [~, dist2] = cellfun(distfun2,fg(jj).fibers,'UniformOutput',false);
    end
    % record which fibers need to be flipped. We flip any fiber who's
    % las node is closer to the starting ROI than its first node
    flipped{jj} = cellfun(@(x) x(1)>x(2), dist1,'UniformOutput',false);
    % Flip fibers that need to be flipped CHECK THIS!!!
    fg(jj).fibers = cellfun(@flipfun,fg(jj).fibers,flipped{jj},'UniformOutput',false);
    % Identify fibers that don't have a startpoint and endpoint within 4mm
    % of the defining ROIs
    if ~isempty(Lnum2)
        keep{jj} = cellfun(@(x1,x2) min(x1) < dCritSq && min(x2) < dCritSq,dist1,dist2);
    else
        keep{jj} = cellfun(@(x1) min(x1) < dCritSq, dist1);
    end
    % Discard the fibers
    fg(jj).fibers = fg(jj).fibers(keep{jj});
end

% Convert fiber group structure to whatever type was input
if array == 0
    fg = dtiFgArrayToFiberGroup(fg,fgName);
end

return

function fibOut = flipfun(fibIn,ind)
% Function to flip fibers
if(ind ==1)
    fibOut = fliplr(fibIn);
else
    fibOut = fibIn;
end

return

function Lnum = getLabelNumber(Lname)

if isempty(Lname)
    Lnum = [];
    return
end
% Get the AFQ base directory
AFQbase = AFQ_directories;
% Template directory
tdir = fullfile(AFQbase,'templates','labelMaps');
% Load template region names
labels = readTab(fullfile(tdir,'MNI_AAL.txt'),',');
% Get the desired label number for the startpoint
Lnum = find(strcmpi(Lname,labels(:,2)));
% If the name was not in the label file check if it is a composite region
if isempty(Lnum)
    switch(Lname)
        case'occipital'
            % We do not include the fusiform
            Lnum = 43:54;
        case 'leftoccipital'
            Lnum = [43:2:53];
        case 'rightoccipital'
            Lnum = [44:2:54];
        case'temporal'
            % We include the fusiform
            Lnum = [37:42 55, 56, 79:90];
        case 'lefttemporal'
            [37:2:41 55 79:2:89];
        case 'righttemporal'
            [38:2:42 56 80:2:90];
        case 'leftanttemporal'
            % 41 = amygdala
            [41 83 87];
        case 'rightanttemporal'
            [42 84 88];
        case 'leftfrontal'
            [1:2:25];
        case 'rightfrontal'
            [2:2:26];
        case 'leftparietal'
            [57:67];
        case 'rightparietal'
            [58:68];
        case 'leftinfparietal'
            [61 63 65];
        case 'rightinfparietal'
            [62 64 65];
        case 'cerebellum'
            91:116;
    end
end

