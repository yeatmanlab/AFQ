function msh = AFQ_meshCreate(im, varargin)
% Create a mesh structure and build a 3d mesh from a nifti image
%
% msh = AFQ_meshCreate(im)
% msh = AFQ_meshCreate(im, 'smooth')

%% Allocate the fields of the mesh structure

% mesh.tr corresponds to the structure expected by matlab's patch function.
msh.type                = 'afq mesh version 1';
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
    
    % Smooth the vertices and add them to the structure. Start with 20
    % smoothing iterations and work our way up. The smooth mesh function
    % will act on the vertices stored in msh.vertices. So after each
    % smoothing iteration we add the new vertices to this field and the
    % perform the operation again
    tmp = smoothpatch(msh.tr, [], 20);
    msh.vertex.smooth20 = tmp.vertices;
    
    % Only continue smoothing if desired
    sm = find(strcmpi('smooth',varargin));
    if ~isempty(sm)
        for ii = [20 40 80 160 320]
            tmp = smoothpatch(tmp, [], ii);
            eval(sprintf('msh.vertex.smooth%d = tmp.vertices;',ii*2))
        end
    end
    
    %% Color the mesh vertices
    
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
        
    elseif isfield(params,'color')
        % Or if a color was defined for the mesh set that as the base color
        msh = AFQ_meshSet(msh,'basecolor', params.color);
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
