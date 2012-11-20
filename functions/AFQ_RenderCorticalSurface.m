function p = AFQ_RenderCorticalSurface(segIm, color, a)
% Render the cortical surface from a binary segmentation image
%
% h = AFQ_RenderCorticalSurface(segIm)
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
%
% Outputs:
% p       - Handel for the patch object that was added to the figure
%           window. The rendering can be deleted with delete(p)
%
% Copyright Jason D. Yeatman November 2012

if ~exist('color','var') || isempty(color)
    color = [.8 .7 .6];
end
if ~exist('a','var') || isempty(alpha)
    a = 1;
end
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
% Render the cortical surface
p = patch(tr,'facecolor',color,'edgecolor','none');
alpha(p,a);
axis('image');axis('vis3d');
set(p,'specularstrength',.5,'diffusestrength',.75);