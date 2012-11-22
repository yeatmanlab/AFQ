function msh = AFQ_meshCreate(im, varargin)
% Create a mesh structure and build a 3d mesh from a nifti image
%
% msh = AFQ_meshCreate(im)
% msh = AFQ_meshCreate(im, 'smooth')

%% Allocate the fields of the mesh structure
msh.vertices         = [];
msh.faces            = [];
msh.FaceVertexCData  = [];
msh.vertex.origin    = [];
msh.vertex.smooth20  = [];
msh.vertex.smooth40  = [];
msh.vertex.smooth80  = [];
msh.vertex.smooth160 = [];
msh.vertex.smooth320 = [];
msh.vertex.smooth640 = [];

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
    
    % Add the faces to the msh structure
    msh.faces = tr.faces;
    
    % Smooth the vertices and add them to the structure. Start with 20
    % smoothing iterations and work our way up. The smooth mesh function
    % will act on the vertices stored in msh.vertices. So after each
    % smoothing iteration we add the new vertices to this field and the
    % perform the operation again
    msh.vertices = msh.vertex.origin;
    tmp = smoothpatch(msh, [], 20);
    msh.vertex.smooth20 = tmp.vertices; 
    
    % Only continue smoothing if desired
    sm = find(strcmpi('smooth',varargin));
    if ~isempty(sm)
        for ii = [20 40 80 160 320]
            tmp = smoothpatch(tmp, [], ii);
            eval('msh.vertex.smooth%d = tmp.vertices;',ii*2)
        end
    end
    
    % Set the default vertices to render as the vertices with 20 smoothing
    % iterations
    msh = AFQ_meshSet(msh, 'vertices', 'smooth20');
    
    end
end