function h = AFQ_RenderRoi(roi, color, method, render, varargin)
% Render an ROI as a 3D surface
%
% h = AFQ_RenderRoi(roi, color , [method = 'trimesh'], [render = 'surface'])
%
% Inputs:
%
% roi    = Roi structure or an Nx3 matrix of X,Y,Z coordinates
% color  = The color to render the roi. Default is red: color = [1 0 0]
% method = There are many ways to build a surface mesh from coordinates.
%          The default is method = 'mesh' which uses the AFQ_meshCreate
%          function to build and AFQ mesh and render it. This is by far the
%          best and other options are only left in for legacy.
%          method = 'isosurface' converts the coordinates to an image than
%          uses isosurface to build a mesh. method = 'trimesh' builds a
%          triangle mesh out of the coordinates and redners that as a
%          surface mesh.
% render = What to render. Either the roi surface: render = 'surface' or a
%          wire frame: render = 'wire' 
%
% Example: fg = dtiReadFibers('L_Arcuate.mat');
% roi = dtiReadRoi('roi1.mat'); AFQ_RenderFiber(fg); % Render the fibers
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
if ~exist('method', 'var') || isempty(method)
    method = 'mesh';
end
if ~exist('render','var') || isempty(render)
    render = 'surf';
end

% % Remove internal coordinates from the ROI so that all that are
% % left are the surface coordinates.
% % Convert coordinates to an image
% [roiImg, imgXform, bb] = dtiRoiToImg(coords);
% % Define the perimeter of the image
% roiImg = bwperim(roiImg);
% % Convert the image back to x,y,z coordinates
% roi = dtiRoiFromImg(roiImg, inv(imgXform), bb);
% coords = roi.coords;

%% Render the ROI

% If an variables are in varargin put them in an afq params structure
params = CreateParamsStruct(varargin);

% Get the current figure window handle
f = gcf;
% Check if there is anything plotted in the figure window already and only
% open a new plotting window if necesary
newfig = isempty(get(f,'children'));

% Compute the range of the coordinates in the roi
roi_min = min(coords);
roi_max = max(coords);

% Choose which method to use for rendering
switch(method)
    case {'mesh'}
        
        % Convert the array of coordinates to an image
        [roiImg, imgXform] = dtiRoiToImg(roi);
        
        % Set ROI as a nifti struct and save
        ni = niftiCreate('data',uint8(roiImg),'qto_xyz',imgXform);
        % ni = niftiCreate(uint8(roiImg),imgXform);%Old code
        
        % Set up a parameters structure
        params.color = color; params.newfig = newfig; params.boxfilter = 3;
        
        % Close the figure window that was opened because we will open a
        % new one with a number of properties set how we like
        if params.newfig==1,close(f);end
        
        % Build a mesh from the nifti image of the roi
        h = AFQ_RenderCorticalSurface(ni, params);
        if strcmp(render,'wire')
            warning('\nWire mesh not implimented yet');
        end  
        
    case {'isosurface', 'isosurf'}
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
    case {'trimesh' 'triangle' 'delaunay' 'trianglemesh'}
        % Build the triangle mesh from the coordinates. Matlab changed the
        % delaunay function at some point so first check whether it is the
        % new or old version of the function
        if exist('delaunay','builtin')
            Tri = delaunay(coords);
        else
            Tri = delaunay3(coords(:,1),coords(:,2),coords(:,3));
        end
        switch(render)
            case {'surface' 'surf'}
                % Render the mesh as a surface
                trisurf(Tri, coords(:,1), coords(:,2), coords(:,3),...
                    'facecolor',color,'edgecolor','none',...
                        'FaceAlpha',1,'FaceLighting','gouraud', ...
                        'SpecularColorReflectance',.2,'ambientstrength',.6);
            case {'wire'}
                % Render a wire mesh
                trimesh(Tri, coords(:,1), coords(:,2), coords(:,3),...
                    'edgecolor',color);
        end
end

% Old axis limits
ax = axis;
% If the axis only has X and Y limits set, then add Z limits
if length(ax) == 4
    ax(5) = roi_min(3); ax(6) = roi_max(3);
end
% New axis scaling
ax2 = [roi_min(1) roi_max(1) roi_min(2) roi_max(2) roi_min(3) roi_max(3)];
% Check if any part of the ROI will not be visible in the current axis and
% rescale the axis if needed.
axnew(1:2:6) = min(vertcat(ax(1:2:6), ax2(1:2:6)));
axnew(2:2:6) = max(vertcat(ax(2:2:6), ax2(2:2:6)));
axis(axnew);

% Turn hold on in case other features are added to the rendering
hold on;