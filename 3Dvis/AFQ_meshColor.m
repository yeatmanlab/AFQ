function msh = AFQ_meshColor(msh, varargin)
% Color the a mesh based on an overlay image or surface curvature
%
% msh = AFQ_meshColor(msh, 'overlay', overlayIm, 'thresh', [min max], 'crange', [min max], 'cmap', [cmap = 'jet'])
%
% AFQ_meshColor takes in a msh mesh structure and adds a new colors field
% to it that specifies an rgb value for each mesh vertex. The colors can be
% computed from an overlay image to heatmap the mesh based on fMRI signal
% or quantitative MRI parameters. Other options are to color the mesh based
% on surface curvature. See help AFQ_RenderCorticalSurface and
% help AFQ_meshCreate for a full description of all the parameters that can
% be set.
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
    if ~isfield(params,'interp') || isempty(params.interp)
        interpMethod = 'spline';
    else
        interpMethod = params.interp;
    end
    cvals = dtiGetValFromImage(overlayIm.data, AFQ_meshGet(msh, 'vertexorigin'), overlayIm.qto_ijk, interpMethod);
    % Remove file extension and path to get the name of the image
    [~,valname] = fileparts(overlayIm.fname);
    % Remove a secondary extension if there is one
    valname = prefix(valname);
    % Set these values to the vals field of the msh structure
    msh = AFQ_meshSet(msh, 'vals', valname, cvals);
    % Find which vertices do not surpass the overlay threshold
    if isfield(params,'thresh') && ~isempty(params.thresh) && length(params.thresh) == 1
        subthresh = cvals < params.thresh;
    elseif isfield(params,'thresh') && ~isempty(params.thresh) && length(params.thresh) == 2
        subthresh = cvals < params.thresh(1) | cvals > params.thresh(2);
    end
    % Convert the values to rgb colors by associating each value with a
    % location on the colormap
    FaceVertexCData = vals2colormap(cvals,cmap,crange);
    % If a threshold was passed in then reasign the default cortex color to
    % vertices that are outside the range defined by threh
    if exist('subthresh','var')
        basecolor = AFQ_meshGet(msh, 'basecolor');
        FaceVertexCData(subthresh,1) = basecolor(1);
        FaceVertexCData(subthresh,2) = basecolor(2);
        FaceVertexCData(subthresh,3) = basecolor(3);
    end
    % Add these colors to the msh structure
    msh = AFQ_meshSet(msh, 'colors', valname, FaceVertexCData);
    % Set this to be the default color for rendering the mesh
    msh = AFQ_meshSet(msh, 'FaceVertexCData', valname);
end

%% Combine multiple color maps if requested
if isfield(params,'combine') && ~isempty(params.combine)
    % Loop over all the requested color maps to combine
    if iscell(params.combine)
        cnames = params.combine;
    else
        cnames{1} = params.combine;
    end
    rgbvals = [];
    for ii = 1:numel(cnames)
        rgbvals(:,:,ii) = msh.colors.(cnames{ii});
        % Find vertices where the colormap has the base color
        base(:,ii) = ismember(rgbvals(:,:,ii),msh.colors.base,'rows');
    end
    
    % Average the colors
    for ii = 1:size(rgbvals,1)
        if sum(~base(ii,:))==0
            cdata(ii,:) = msh.colors.base;
        else
            cdata(ii,:) = mean(rgbvals(ii,:,~base(ii,:)),3);
        end
    end
    % Add this new colormap to the mesh structure
    msh = AFQ_meshSet(msh, 'colors', horzcat(cnames{:}), cdata);
    % Set this to be the default color for rendering the mesh
    msh = AFQ_meshSet(msh, 'FaceVertexCData', horzcat(cnames{:}));
end

return
