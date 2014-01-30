function msh = AFQ_meshCreate(im, varargin)
% Create a mesh structure and build a 3d mesh from a nifti image
%
% msh = AFQ_meshCreate(im)
%
% AFQ_meshCreat will create the AFQ msh mesh structure from a nifti image.
% The variable im can either be a path to a nifti image or a nifti image
% file. If no image is passed in then an empty strucure will be returned.
% The mesh can then be rendered with AFQ_RenderCorticalSurface. There are a
% number of parameters that can be set with the form of:
% msh = AFQ_meshCreate(im, 'param', value). 
% All the examples should run if you first run:
% [~, AFQdata] = AFQ_directories; cd(fullfile(AFQdata,'mesh'));
% im = 'segmentation.nii.gz'; overlayIm = 'Left_Arcuate_Endpoints.nii.gz';
%
% Examples:
%
% % Create a mesh structure where the cortex is colored a solid color (other
% % than "brain" color)
% msh = AFQ_meshCreate(im, 'color', [color = [.8 .7 .6]])
%
% % Create a mesh structure with a set of vertices that are constructed from 
% % a filtered version of the image ('box' is the only implimented filter).
% % This will create a cortical surface that looks fuller -- like the gray 
% % matter is thicker 
% msh = AFQ_meshCreate(im ,'boxfilter', 5)
% 
% % Create a mesh structure and compute vertices with various numbers of
% % smoothing iterations. More smoothing iterations makes a more inflated,
% % less gyrified cortical surface. Set of vertices with all of the listed
% % numbers of smoothing iterations will be saved in the mesh structure. This
% % way you can render the mesh with a different amount of smoothing by
% % simply setting those to be the default vertices with a call like:
% % AFQ_meshSet(msh, 'vertices', 200)
% msh = AFQ_meshCreate(im, 'smooth', [20 40 80 200])
% 
% % Create a mesh structure with a set of mesh vertice colors computed from
% % an overlay image. The coloring can be computed with various thresholds
% % ('thresh') meaning that values below the defined min and above the
% % defined max will be colored "cortex" color and not heatmapped based on
% % the overlay image. The color range ('crange') can also be set if you want
% % to have the heatmap saturate below the maximum value or above the minimum
% % value. You can also change the colormap ('cmap'). 
% % The default is jet (see doc colormap).
% msh = AFQ_meshCreate(im, 'overlay', overlayImage, 'thresh', [min max],...
% 'crange', [min max], 'cmap', [cmap = 'jet'])
%
% Copyright Jason D. Yeatman November 2012

%% Allocate the fields of the mesh structure

% mesh.tr corresponds to the structure expected by matlab's patch function.
msh.type                = 'afq mesh version 1';
msh.image               = [];
msh.tr.vertices         = [];
msh.tr.faces            = [];
msh.tr.FaceVertexCData  = [.8 .7 .6];
% msh.vertex saves the various vertices with different smoothing
% parameters
msh.vertex.current   = 'origin';
msh.vertex.origin    = [];
msh.vertex.smooth20  = [];
% msh.face saves the faces corresponding to each set of vertices. These are
% generally the same
msh.face.current     = 'origin';
msh.face.origin      = [];
msh.face.smooth20    = 'origin';
% Values associated with mesh vertices
msh.vals             = [];
% Colors associated with mesh vertices
msh.colors.current   = 'base';
msh.colors.base      = [.8 .7 .6]; % Base color of mesh
% A field for ROIs associated with mesh vertices
msh.rois             = [];
msh.roi.show         = {};
% A field for fibers associated with mesh vertices
msh.fibers           = [];
msh.fibers.show      = {};
% Field to store affine and other transormations of mesh vertices
msh.xform            = [];
%% Colect parameters from varargin and put them in a structure
if length(varargin) == 1 && isparams(varargin{1})
    params = varargin{1};
else
    params = CreateParamsStruct(varargin);
end

%% Build a mesh if an image was passed in
if exist('im','var') && ~isempty(im)
    
    % Load the image
    if ischar(im), im = readFileNifti(im); end
    
    % Save the image path in the mesh structure
    msh.image = im.fname;
    
    % Add the image affines to the mesh structure
    msh.xform.qto_ijk = im.qto_ijk;
    msh.xform.qto_xyz = im.qto_xyz;
    
    % permute the image dimensions (This is because the x,y and z dimensions in
    % matlab do not correspond to left-right, anterior-posterior, up-down.
    data = permute(im.data, [2 1 3]);
    
    % Check if the mesh is in mrvista format. If so convert the white
    % matter in the mesh into a binary image
    if sum(data(:)==3)>0 && sum(data(:)==4)>0
        data = uint8(data==3 | data==4);
    end
    
    % Build the mesh
    tr = isosurface(data,.1);
    
    % Transform the vertices of the mesh to acpc space. 
    tr.vertices = mrAnatXformCoords(im.qto_xyz,tr.vertices);
    
    % Reduce the number of faces for computational efficiency
     if isfield(params, 'reduce') && ~isempty(params.reduce)
        tr = reducepatch(tr, params.reduce);
     end
    % Add these vertices and faces to the mesh structure. These are the
    % original vertices and will be used in the future to referece
    % coordinates of overlay images
    msh = AFQ_meshSet(msh,'vertex','origin',tr.vertices,tr.faces);
    msh = AFQ_meshSet(msh, 'vertices','origin');
    
    %% Smooth the mesh
    
    % Smooth the vertices and add them to the structure. Start with 20
    % smoothing iterations and work our way up. The smooth mesh function
    % will act on the vertices stored in msh.vertices. So after each
    % smoothing iteration we add the new vertices to this field and the
    % perform the operation again.
    if ~isfield(params,'smooth') || isempty(params.smooth)
        % 20 smoothing iterations is the default
        params.smooth = 20;
    end
    numiter = 0;
    % Get the triangle mesh so we can start smoothing it
    tmp = AFQ_meshGet(msh,'triangles');
    for ii = params.smooth
        % Keep track of the number of iterations that have already been
        % done so we don't have to recompute what was already finished
        numiter = ii - numiter;
        % Perform the smoothing
        tmp = smoothpatch(tmp, [], numiter);
        % Add these data to the mesh structure in a properly named field
        name = sprintf('smooth%d', ii);
        msh = AFQ_meshSet(msh, 'vertex', name, tmp.vertices);
    end
    
    %% Create a mesh from a smoothed version of the segmentation image
    if isfield(params, 'boxfilter') && ~isempty(params.boxfilter)
        fname = sprintf('filteredbox%d',params.boxfilter);
        % smooth the image with the desired filter size
        data = smooth3(data,'box',params.boxfilter);
        % make a mesh
        tr = isosurface(data,.1);
        % transform the vertices to acpc space
        tr.vertices = mrAnatXformCoords(im.qto_xyz, tr.vertices);
        % smooth the vertices with 20 smoothing iterations
        tr = smoothpatch(tr,[],20);
        % Find the correspoding vertices in the unfiltered image. This is
        % essential for propperly mapping data to the surface of this mesh
        tr.map2origin = nearpoints(tr.vertices', msh.vertex.origin');
        % Save this within the mesh structure
        msh = AFQ_meshSet(msh, 'vertex',fname, tr.vertices, tr.faces, tr.map2origin); 
        % Set these to be the vertices that are rendered
        msh = AFQ_meshSet(msh, 'vertices', fname);
    end
    
    %% Color the mesh vertices
    
    % Set the base color of the mesh if it was defined
    if isfield(params,'color')
        % Or if a color was defined for the mesh set that as the base color
        msh = AFQ_meshSet(msh,'basecolor', params.color);
    end
    
    % If an overlay image was provided then use that to color the mesh
    if isfield(params,'overlay') && ~isempty(params.overlay)
        msh = AFQ_meshColor(msh, params);
    end
    
    %% Set the default vertices and color for mesh rendering
    
    % By default we render with 20 smoothing iterations
    if ~isfield(params, 'boxfilter') || isempty(params.boxfilter)
        msh = AFQ_meshSet(msh, 'vertices', ['smooth' num2str(params.smooth)]);
    end
    % If no overlay was sent in then color each vertex the base color
    if ~isfield(params,'overlay') || isempty(params.overlay) 
        msh = AFQ_meshSet(msh, 'FaceVertexCData', 'base');
    end
    
end
