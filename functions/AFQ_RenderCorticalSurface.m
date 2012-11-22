function p = AFQ_RenderCorticalSurface(segIm, color, a, overlay, thresh, crange, cmap, newfig)
% Render the cortical surface from a binary segmentation image
%
% p = AFQ_RenderCorticalSurface(segIm, color, a, overlay, thresh, crange, cmap, newfig)
%
% This function takes in a segmentation image and renders it in 3d. It is
% optimized to look good for the cortical surface but any image will work.
% The rendering will be added to the current figure window so you can add
% the cortex to a rendering of fiber groups and adjust it's transparency
%
% Inputs:
% segIm   - A path to a nifti image to render. It must be a binary mask
% color   - RGB value for the surface of the rendering. Default is "brain"
%           color
% a       - The transparency of the surface (alpha). 0 is completely
%           transparent and 1 is completely opaque
% overlay - Another image to use to color the surface (for example an fMRI
%           contrast).
% thresh  - A threshold above/below which now overlay values will be
%           painted on the cortex and the cortex will be left cortex
%           colored. Thresh can be a single number (minumum) or a vector of
%           2 numbers (minimum and maximum).
% crange  - Define which overlay values should be mapped to the minimum and
%           maximum values of the color map. All values below crange(1)
%           will be colored the minimum value and all values above
%           crange(2) will be colored the maximum value. The default color
%           range is defined by the minimum and maximum values of the
%           overlay image that get mapped to any mesh vertex.
% cmap    - Name of the colormap to use. For example 'jet' or 'autumn'
% newfig  - Whether or not to open a new figure window
%
% Outputs:
% p       - Handel for the patch object that was added to the figure
%           window. The rendering can be deleted with delete(p)
%
% Example:
%
% % Get data
% [~, AFQdata] = AFQ_directories; 
% segIm = fullfile(AFQdata,'mesh','segmentation.nii.gz');
% overlay = fullfile(AFQdata,'mesh','t1.nii.gz');
% % Render the cortical surface colored by the T1 values at each vertex
% p = AFQ_RenderCorticalSurface(segIm, [], [], overlay)
%
% Copyright Jason D. Yeatman November 2012

if ~exist('color','var') || isempty(color)
    color = [.8 .7 .6];
end
if ~exist('a','var') || isempty(a)
    a = 1;
end
% Read the overlay image if a path was provided
if exist('overlay','var') && ~isempty(overlay) && ischar(overlay)
    overlay = readFileNifti(overlay);
end   
% Default color map is jet
if ~exist('cmap','var') || isempty(cmap)
    cmap = 'jet';
end
% Default color range is defined by the values in the overlay
if ~exist('crange','var') || isempty(crange)
    crange = [];
end
if ~exist('newfig','var') || isempty(newfig)
    newfig = 1;
end
%% Build a mesh of the cortical surface

% Load the image
im = readFileNifti(segIm);
% permute the image dimensions (This is because the x,y and z dimensions in
% matlab do not correspond to left-right, anterior-posterior, up-down.
data = permute(im.data, [2 1 3]);
% smooth the image
data = smooth3(data,'box',5);
% make a mesh
tr = isosurface(data,.1);
% transform the vertices to acpc space
tr.vertices = mrAnatXformCoords(im.qto_xyz,tr.vertices);

%% Color the mesh vertices

% If an overlay image was provided use that to color the mesh, otherwise
% color it all a uniform color
if exist('overlay','var') && ~isempty(overlay)
    % Interpolate overlay values at each vertex of the mesh
    cvals = dtiGetValFromImage(overlay.data, tr.vertices, overlay.qto_ijk, 'spline');
    % Find which vertices do not surpass the overlay threshold
    if exist('thresh','var') && ~isempty(thresh) && length(thresh) == 1
        subthresh = FaceVertexCData < thresh;
    elseif exist('thresh','var') && ~isempty(thresh) && length(thresh) == 2
        subthresh = FaceVertexCData < thresh(1) || FaceVertexCData > thresh(2);
    end
    % Convert the values to rgb colors by associating each value with a
    % location on the colormap
    tr.FaceVertexCData = vals2colormap(cvals,cmap,crange);
    % If a threshold was passed in then reasign the default cortex color to
    % vertices that are outside the range defined by threh
    if exist('subthresh','var')
        tr.FaceVertexCData(subthresh,:) = color;
    end
else
    tr.FaceVertexCData = repmat(color,size(tr.vertices,1),1);
end
%% Render the cortical surface
if newfig == 1
    figure;
end
% Use patch to render the mesh
p = patch(tr);
%p = patch(tr,'facecolor',color,'edgecolor','none');

% Interpolate the coloring along the surface
shading('interp'); 
% Set the type of lighting
lighting('gouraud');
% Set the alpha
alpha(p,a);
% Set axis size
axis('image');axis('vis3d');
% Set lighiting options of the cortex
set(p,'specularstrength',.5,'diffusestrength',.75);

% If it was a new figure window add a light to it
if newfig == 1
    camlight('right');
end
