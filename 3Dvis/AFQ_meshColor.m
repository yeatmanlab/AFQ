function msh = AFQ_meshColor(msh, varargin)
% Color the a mesh based on an overlay image or surface curvature
%
% msh = AFQ_meshColor(msh, 'overlay', overlayIm, 'thresh', [min max], 'crange', [min max])
%
% AFQ_meshColor takes in a msh mesh structure and adds a new colors field
% to it that specifies an rgb value for each mesh vertex. The colors can be
% computed from an overlay image to heatmap the mesh based on fMRI signal
% or quantitative MRI parameters. Other options are to color the mesh based
% on surface curvature.
%
% Arguments can be passed in in one of 2 ways. (1) In the form:
%'parameter', [value] - where the parameter is defined and then the
% arguments or inputs follow the parameter name. (2) As a parameters
% structure - see CreateParamsStruct.

% Put all the arguments into a params structure
if length(varargin) == 1 && isparams(varargin{1})
    params = varargin{1};
else
    params = CreateParamsStruct(varargin);
end

%% If an overlay image was provided use that to color the mesh
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
    if isfield(params,'thresh') && ~isempty(params.thresh) && length(params.thresh) == 1
        subthresh = cvals < params.thresh;
    elseif isfield(params,'thresh') && ~isempty(params.thresh) && length(params.thresh) == 2
        subthresh = cvals < params.thresh(1) || cvals > params.thresh(2);
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
    % Set this to be the default color for rendering the mesh
    msh = AFQ_meshSet(msh, 'FaceVertexCData', valname);
end