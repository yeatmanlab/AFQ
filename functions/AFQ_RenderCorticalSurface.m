function p = AFQ_RenderCorticalSurface(segIm, color, a)
% Render the cortical surface from a binary segmentation image
%
% h = AFQ_RenderCorticalSurface(segIm)

if ~exist('color','var') || isempty(color)
    color = [.8 .7 .6];
end
if ~exist('a','var') || isempty(alpha)
    a = 1;
end
% Load the image
im = readFileNifti(segIm);
% permute the image dimensions
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
set(p,'specularstrength',.5,'diffusestrength',.75)