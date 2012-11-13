function [fg keep flipped] = AFQ_DefineFgEndpoints(fg,startpoint,endpoint,dt6Path,dCrit,invDef)
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
    dCritSq = 16;
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
% Number of fiber groups
nFG = length(fg);
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
    if iscell(startpoint) && iscell(endpoint)
        Lnum1 = getLabelNumber(startpoint{jj});
        Lnum2 = getLabelNumber(endpoint{jj});
    else
        Lnum1 = getLabelNumber(startpoint);
        Lnum2 = getLabelNumber(endpoint);
    end
    % Make ROIs in native space from the template ROIs Load the desired
    % startpoint ROI from the AAL template and transform it
    % into the individual's native space
    [~, invDef, roiStart] = dtiCreateRoiFromMniNifti(dt6Path,AAL,invDef,0,Lnum1);
    [~,~, roiEnd]  = dtiCreateRoiFromMniNifti(dt6Path,AAL,invDef,0,Lnum2);
    % Make a function to compute the distance between fiber endpoints and
    % the desired start and endpoint
    distfun1 = @(x) nearpoints([x(:,1) x(:,end)],roiStart.coords');
    distfun2 = @(x) nearpoints([x(:,1) x(:,end)],roiEnd.coords');
    % Compute distance for each fiber to the start and end ROIs
    [~, dist1] = cellfun(distfun1,fg(jj).fibers,'UniformOutput',false);
    [~, dist2] = cellfun(distfun2,fg(jj).fibers,'UniformOutput',false);
    % record which fibers need to be flipped
    flipped{jj} = cellfun(@(x) x(1)>x(2), dist1,'UniformOutput',false);
    % Flip fibers that need to be flipped CHECK THIS!!!
    fg(jj).fibers = cellfun(@flipfun,fg(jj).fibers,flipped{jj},'UniformOutput',false);
    % Discard fibers that don't have a startpoint and endpoint within 4mm
    % of the defining ROIs
    keep{jj} = cellfun(@(x1,x2) min(x1) < dCritSq && min(x2) < dCritSq,dist1,dist2);
    fg(jj).fibers = fg(jj).fibers(keep{jj});
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
        case('occipital')
            % We include the fusiform
            Lnum = 43:56;
        case{'temporal'}
            % We include the fusiform
            Lnum = [55, 56, 79:90];
    end
end

