function msh = AFQ_meshCreate(im, varargin)
% Create a mesh structure and build a 3d mesh from a nifti image
%
% msh = AFQ_meshCreate(im)
% msh = AFQ_meshCreate(im ,'boxfilter', 5)
% msh = AFQ_meshCreate(im, 'smooth', [20 40 80 200])
% msh = AFQ_meshCreate(im, 'overlay', overlayImage)
% msh = AFQ_meshCreate(im, 'color', [.8 .7 .6])

%% Allocate the fields of the mesh structure

% mesh.tr corresponds to the structure expected by matlab's patch function.
msh.type                = 'afq mesh version 1';
msh.image               = [];
msh.tr.vertices         = [];
msh.tr.faces            = [];
msh.tr.FaceVertexCData  = [.8 .7 .6];
% mesh.vertex saves the various vertices with different smoothing
% parameters
msh.vertex.origin    = [];
msh.vertex.smooth20  = [];
% Values associated with mesh vertices
msh.vals             = [];
% Colors associated with mesh vertices
msh.colors.base      = [.8 .7 .6]; % Base color of mesh

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
    
    % permute the image dimensions (This is because the x,y and z dimensions in
    % matlab do not correspond to left-right, anterior-posterior, up-down.
    data = permute(im.data, [2 1 3]);
    
    % Build the mesh
    tr = isosurface(data,.1);
    
    % Transform the vertices of the mesh to acpc space and add them to our
    % msh structure. These are the original vertices and will be used in
    % the future to referece coordinates of overlay images
    msh.vertex.origin = mrAnatXformCoords(im.qto_xyz,tr.vertices);
    msh.tr.vertices = msh.vertex.origin;
    
    % Add the faces to the msh structure
    msh.tr.faces = tr.faces;
    
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
    for ii = params.smooth
        % keep track of the number of iterations that have already been
        % done so we don't have to recompute what was already finished
        numiter = ii - numiter;
        tmp = smoothpatch(tmp, [], numiter);
        eval(sprintf('msh.vertex.smooth%d = tmp.vertices;',ii))
    end
    
    %% Create a mesh from a smoothed version of the segmentation image
    if isfield(params, 'boxfilter') || ~iesmpty(params.boxfilter)
        fname = sprintf('box%d',params.boxfilter);
        % smooth the image with the desired filter size
        data = smooth3(data,'box',params.boxfilter);
        % make a mesh
        tr = isosurface(data,.1);
        % transform the vertices to acpc space
        tr.vertices = mrAnatXformCoords(im.qto_xyz,msh.vertices);
        % Find the correspoding vertices in the unfiltered image. This is
        % essential for propperly mapping data to the surface of this mesh
        tr.map2origin = nearpoints(tr.vertices, msh.vertex.origin);
        % Save this within the mesh structure
        AFQ_meshSet(msh, 'filter',fname, tr)     
    end
    
    %% Color the mesh vertices
    
    % Set the base color of the mesh if it was defined
    if isfield(params,'color')
        % Or if a color was defined for the mesh set that as the base color
        msh = AFQ_meshSet(msh,'basecolor', params.color);
    end
    
    % If an overlay image was provided use that to color the mesh, otherwise
    % color it all a uniform color
    if isfield(params,'overlay') && ~isempty(params.overlay)
        % Load the image if overlay is a path
        if ischar(params.overlay)
            overlayIm = readFileNifti(params.overlay);
        else
            overlayIm = params.overlay;
        end
        % Default color map is jet
        if ~isfield(params,'cmap')
            cmap = 'jet';
        else
            cmap = params.cmap;
        end
        % Default color range is defined by the values in the overlay
        if ~isfield(params,'crange')
            crange = [];
        else 
            crange = params.crange;
        end
        % Interpolate overlay values at each vertex of the mesh
        cvals = dtiGetValFromImage(overlayIm.data, AFQ_meshGet(msh, 'vertexorigin'), overlayIm.qto_ijk, 'spline');
        % Remove file extensions to get the name of the image
        valname = prefix(prefix(overlayIm.fname));
        % Set these values to the vals field of the msh structure
        msh = AFQ_meshSet(msh, 'vals', valname, cvals);
        % Find which vertices do not surpass the overlay threshold
        if exist('thresh','var') && ~isempty(thresh) && length(thresh) == 1
            subthresh = cvals < thresh;
        elseif exist('thresh','var') && ~isempty(thresh) && length(thresh) == 2
            subthresh = cvals < thresh(1) || cvals > thresh(2);
        end
        % Convert the values to rgb colors by associating each value with a
        % location on the colormap
        FaceVertexCData = vals2colormap(cvals,cmap,crange);
        % If a threshold was passed in then reasign the default cortex color to
        % vertices that are outside the range defined by threh
        if exist('subthresh','var')
            FaceVertexCData(subthresh,:) = AFQ_meshGet(msh, 'basecolor');
        end
        % Add these colors to the msh structure
        msh = AFQ_meshSet(msh, 'colors', valname, FaceVertexCData);
    end
    
    %% Set the default vertices and color for mesh rendering
    
    % By default we render with 20 smoothing iterations
    msh = AFQ_meshSet(msh, 'vertices', 'smooth20');
    % Set the color to be that of the overlay
    if exist('overlayIm','var') && ~isempty(overlayIm)
        msh = AFQ_meshSet(msh, 'FaceVertexCData', valname);
    else
        % If no overlay was sent in the set the default color to every
        % vertex
        msh = AFQ_meshSet(msh, 'FaceVertexCData', 'base');
    end
    
end
