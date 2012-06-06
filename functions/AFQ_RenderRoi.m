function h = AFQ_RenderRoi(roi, color)
% Render an ROI as a 3D surface
%
% h = AFQ_RenderRoi(roi, color)
%
% Inputs:
%
% roi    = Roi structure or an Nx3 matrix of X,Y,Z coordinates
% color  = The color to render the roi. Default is red: color = [1 0 0]
%
% Example:
% fg = dtiReadFibers('L_Arcuate.mat'); roi = dtiReadRoi('roi1.mat');
% AFQ_RenderFiber(fg); % Render the fibers
% AFQ_RenderRoi(roi, [1 .5 0]); % Render the roi in orange
%
% Copyright Jason D Yeatman June 2012

%% Check arguments
if ~exist('roi','var') || isempty(roi)
    error('You must supply an ROI')
elseif isstruct(roi)
    coords = roi.coords;
end
if ~exist('color','var') || isempty(color)
    color = [1 0 0];
end

%% Render the ROI
% Compute the range of the coordinates in the roi
roi_min = min(coords);
roi_max = max(coords);
% Create a binary image from the ROI
im = CoordsToImg(coords);
% Specify the corresponding X,Y,Z coordinate of each image index
[X Y Z] = meshgrid(roi_min(1)-1:roi_max(1)+1,roi_min(2)-1:roi_max(2)+1,roi_min(3)-1:roi_max(3)+1);
% Meshgrid outputs the dimensions in a different order so they must be
% permuted
X = permute(X,[2 1 3]);Y = permute(Y,[2 1 3]);Z = permute(Z,[2 1 3]);
% Build a surface mesh of the image
m = isosurface(X,Y,Z,im);
% Render the surface.
h = patch(m,'FaceColor',color,'EdgeColor','none');
% h = isonormals(im,h);