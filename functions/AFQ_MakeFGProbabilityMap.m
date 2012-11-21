function [fgPmap, fgEndPmap] = AFQ_MakeFGProbabilityMap(dtPaths,fgPaths,writeImage)
% Make a fiber tract probability map from fiber group segmentations
%
% [fgPmap, fgEndPmap] = AFQ_MakeFGProbabilityMap(dtPaths, fgPaths, [writeImage = false])
%
% Make a nifti image of a fiber tract probability map from fiber group
% segmentations in a sample of subjects. Each voxel in the resulting
% probability map will represent the number of subjects with fibers in that
% voxel. The resulting nifti images will be in MNI space.
%
% Inputs:
%
% dtPaths   = A cell array of paths to each subject's dt6.mat file or and
%             afq structure. If an afq structure is passed in then paths
%             will be extracted from the structure
% fgPaths   = A cell array of paths to each subject's fiber group or if an
%             afq structure was passed in then simply the name of the fiber
%             group is sufficient.
% writeImage= Give the path to where you would like the probability map to
%             be written (if you would like it to be saved)
%
% Outputs:
%
% fgPmap    = A nifti image of the fiber tract probability map
% fgEndPmap = A nifti image of the fiber tract endpoint probability map
% 
%
% Copyright Jason D. Yeatman, November 2012

if ~exist('fgPaths','var') || isempty(fgPaths)
    error('Please provide paths to the fiber groups')
elseif ischar(fgPaths)
    if isafq(dtPaths)
        % Get the fiber group names from the afq structure and make them
        % lower case
        fgnames = AFQ_get(dtPaths,'fgnames');
        % Format the name of the fiber group to remove spaces etc
        name = mrvParamFormat(fgPaths);
        for jj = 1:length(fgnames)
            fgnames{jj} = mrvParamFormat(dtPaths.fgnames{jj});
        end
        % Check if the desired fiber group is in the structure and get it's
        % number if it is
        fgnum = find(strcmpi(name,fgnames));
    else
       error('Please provide an afq structure or paths to the fiber groups') 
    end
end

% Check the dtPaths. If an afq structure was passed in the get the paths
% out of it
if ~exist('dtPaths','var') || isempty(dtPaths)
    error('please provide paths to each subjects dt6 file');
elseif isafq(dtPaths)
   dtPaths = AFQ_get(dtPaths, 'dt6path'); 
end
if ~exist('writeImage','var') || isempty(writeImage)
    writeImage = 0;
elseif writeImage == 1
    % If no directory is specified save the image in the AFQ data directory
    [~, writeImage] = AFQ_directories;
    fprintf('\nWriting probability map to %s',writeImage)
end

% If the fiber group was one contained in an afq structure that was passed
% in then build the correct fiber group paths
if exist('fgnum','var') && ~isempty(fgnum)
    fgPaths = cell(size(dtPaths));
    for ii = 1:length(dtPaths)
        fgPaths{ii} = fullfile(fileparts(dtPaths{ii}),'fibers', 'MoriGroups.mat');
    end
    % If the fiber group was not in the afq structure then assume it is the
    % name of another fiber group within each subject's fibers directory
elseif exist('fgnum','var') && isempty(fgnum)
    name = fgPaths;
    fgnum = 1;
    for ii = 1:length(dtPaths)
        fgPaths{ii} = fullfile(fileparts(dtPaths{ii}),'fibers', name);
    end
else
    fgnum = 1;
end
% Initialize spm defualts for normalization
spm_defaults; global defaults; params = defaults.normalise.estimate;
% spm_get_defaults - For SPM8
% Set the directory where templates can be found
tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');
% Spatially normalize diffusion data with the MNI (ICBM) template
template = fullfile(tdir,'MNI_JHU_T2.nii.gz');
% Load the template image
tIm = readFileNifti(template);

% Allocate a 4 dimensional array where the 4th dimension corresponds to
% an image of each subject's fiber group
fgIm = zeros([size(tIm.data) length(dtPaths)]);
fgEndIm = zeros([size(tIm.data) length(dtPaths)]);

% loop over each subject
for ii = 1:length(dtPaths)
    %% Load fiber group and dt6
    
    fg = dtiReadFibers(fgPaths{ii});
    % In case fg contains many groups, extract the correct one
    fg = fg(fgnum);
    dt = dtiLoadDt6(dtPaths{ii});
    
    %% Compute spatial normalization
    
    % Rescale b0 image values to get better gary/white/CSF contrast
    alignIm = mrAnatHistogramClip(double(dt.b0),0.3,0.99);
    % Compute normalization
    [sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(alignIm, dt.xformToAcpc, template, params);
    
    %% Transform fibers to MNI space and make an image of the fiber group
    
    % Spatially normalize fibers
    fg_sn = dtiXformFiberCoords(fg, invDef);
    % Concatinate all the fiber coordinates
    fc = horzcat(fg_sn.fibers{:})';
    % Construct another array with the start and endpoints of each fiber
    % group
    fcEnd = [];
    for jj=1:length(fg_sn.fibers)
        fcEnd = vertcat(fcEnd,[fg_sn.fibers{jj}(:,1) fg_sn.fibers{jj}(:,end)]');
    end
    % Remove any coordinates that are nans. This can happen if the
    % coordinates are outside of the MNI image
    fc = fc(~isnan(fc(:,1)),:);
    fcEnd = fcEnd(~isnan(fcEnd(:,1)),:);
    
    % Transform the fiber coordinates to image indices by applying the
    % affine stored in the header of the MNI template image used for
    % normalization. This will shift, scale and rotate the coordinates if
    % necessary
    fc    = mrAnatXformCoords(tIm.qto_ijk, fc);
    fcEnd = mrAnatXformCoords(tIm.qto_ijk, fcEnd);
    % Determine the unique coordinates
    fcU    = unique(round(fc),'rows');
    fcEndU = unique(round(fcEnd),'rows');
    % Convert coordinates to image indices with respect to the template
    % image that was used for normalization
    ind    = sub2ind(size(tIm.data),fcU(:,1),fcU(:,2),fcU(:,3));
    indEnd = sub2ind(size(tIm.data),fcEndU(:,1),fcEndU(:,2),fcEndU(:,3));
    % Create an image of zeros the size of the template image
    tmpIm    = zeros(size(tIm.data));
    tmpEndIm = zeros(size(tIm.data));
    % Put ones at each location where there is a fiber for this subject
    tmpIm(ind)       = 1;
    tmpEndIm(indEnd) = 1;
    % Add this image to the 4th dimension of what will become our fiber
    % probability map
    fgIm(:,:,:,ii)    = tmpIm;
    fgEndIm(:,:,:,ii) = tmpEndIm;
end

%% Create a probability map

% To make sure that we get the transforms correct we will steal all the
% header information from the MNI template image
fgPmap    = tIm;
fgEndPmap = tIm;
% Name it based on the last fiber group that was loaded (but remove capital
% letters and spaces)
fgPmap.fname = mrvParamFormat([fg.name '_probmap.nii.gz']);
fgEndPmap.fname = mrvParamFormat([fg.name '_endpoint_probmap.nii.gz']);

% Count how many subjects had fibers in each voxel by summing across the
% 4th dimension of fgIm. This is our probability map. We replace the data
% within the template image with this new template data.
fgPmap.data    = sum(fgIm,4);
fgEndPmap.data = sum(fgEndIm,4);
% Some visualization programs use this scaling.
fgPmap.cal_max    = max(fgPmap.data(:));
fgEndPmap.cal_max = max(fgEndPmap.data(:));


%% Save the image if desired

if ischar(writeImage)
    curdir = pwd;
    cd(writeImage);
    writeFileNifti(fgPmap);
    writeFileNifti(fgEndPmap);
    cd(curdir);
end
