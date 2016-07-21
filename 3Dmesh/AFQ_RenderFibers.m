function [lightH, fiberMesh, h] = AFQ_RenderFibers(fg,varargin)
% Render a fiber group and tract profile in 3-D
%
%   [lightH. fiberMesh, h] = AFQ_RenderFibers(fg,'PropertyName',PropertyValue, ...)
%
% Given an mrDiffusion fiber group structure, AFQ_RenderFibers(fg), will
% make a 3-D rendering of that fiber group.  The 3-d rendering can be
% rotated with the keyboard arrow keys There are aditional optional
% arguments that can be passed into the function.  They are always defined
% by first entering the 'PropertyName' as a string, followed by the
% property value which may either be a string, a scaler or a data
% structure.  Below all the options are described.
%
% Options:
%
% AFQ_RenderFibers(fg,'dt',dtStruct) - If a dt6 structure is provided
% then the tract FA profile will be rendered along with the fiber group.
% See dtiLoadDt6 and dtiComputeDiffusionPropertiesAlongFG.
%
% AFQ_RenderFibers(fg,'dt',dtStruct,'rois',roi1,roi2) - If a dt6 structure
% and 2 rois are provided then the tract FA profile will onlyh be rendered
% for the portion of the tract spanning between the 2 rois.  See
% dtiReadRois
%
% AFQ_RenderFibers(fg,'dt',dtStruct,'crange',range) - 'crange' defines
% the color range for the tract profile.  The default is range = [0.3 0.6].
% This means that values below 0.3 will take on the minimum value of the
% color map and values above 0.6 will take on the maximum value.
%
% AFQ_RenderFibers(fg,'dt',dtStruct,'cmap',colormap) - 'cmap' defines the
% color map to be used for the tract profile.  The default is colormap =
% 'jet'. Other options are 'autumn' 'hot' 'cool' 'winter' etc. See matlab
% help for other colormaps
%
% AFQ_RenderFibers(fg,'dt',dtStruct,'radius',r) - radius is a 1 by 2 vector
% that defines the radius of the tubes for the fibers and the tract
% profile. r(1) is the radius of the fibers and r(2) is the radius of the
% tract profile. Default is r = [1 5] meaning that fibers have a 1mm radius
% and the tract profile has a 5mm radius.
%
% AFQ_RenderFibers(fg,'tractprofile',TractProfile,'val',['fa']) - A
% tract profile can be passed in if it has been precomputed for the fiber
% tract. TractProfile is a structure created by AFQ_CreateTractProfile and
% contains the coordinates and values to be plotted. The user can also
% define which value in the tract profile should be plotted. The dafault is
% 'fa' but the name of any value that exists within the TractProfile
% structure is ok.
%
% AFQ_RenderFibers(fg,'color', rgbValues) - Render the fiber group in a
% specific color.  rgbValues can be defined in 3 ways to do (1) uniform
% coloring for the fiber group, (2) fiberwise coloring or (3) pointwise
% coloring.
% (1) If rgbValues is a 1 x 3 vector of rgb values than each fiber is
% colored that same color. The default color is gray [0.7 0.7 0.7]. To do
% cyan for example rgbValues = [0 1 1].
% (2) If rgbValues is a N x 3 vector where N is the number of fibers in
% the fiber group, then each fiber in the group is rendered in its own
% color. Each fiber's color is defined by the coresponding row of
% rgbValues. For example to color each fiber a random color:
% rgbValues=rand(length(fg.fibers),3)
% (3) If rgbValues is a 1 x N cell array where N is the number of fibers
% in the fiber group, then each node on fiber n is colored based on the
% corresponding row of rgbValues{n}. This means that each cell must have
% the same number of rows as the corresponding fiber in the fiber group.
% For example to color each point on each fiber based on its FA value:
% vals = dtiGetValFromFibers(dt.dt6,fg,inv(dt.xformToAcpc),'fa');
% rgb = vals2colormap(vals);
% AFQ_RenderFibers(fg,'color',rgb);
% To color each point on each fiber based on values from any nifti image:
% im = readFileNifti('pathToImage.nii.gz');
% vals = dtiGetValFromFibers(im.data,fg,im.qto_ijk)
% rgb = vals2colormap(vals);
% AFQ_RenderFibers(fg,'color',rgb)
%
% AFQ_RenderFibers(fg,'camera', view) - Render the fiber group and view
% from a specific plane.  The options for view are 'sagittal', 'coronal' or
% 'axial'. The default is sagittal. View can also be defined with a 1 by 2
% vector denoting the azimuth (horizontal rotation) and elevation (vertical
% rotation) in degrees of the camera with respect to the fibers. See the
% matlab view function
%
% AFQ_RenderFibers(fg,'jittershading', jitter) - The darkness of each fiber
% can be randomized slightly to make the fiber group take on more of a 3d
% apearance. jitter defines the standard deviation of how amount of
% randomization in fiber darkness.  The default is jitter = 0.2
%
% AFQ_RenderFibers(fg,'jittercolor', jitter) - The color of each fiber can
% be randomized slightly to make each individual streamline stand out more.
% jitter defines the amount of randomization for each color channel with
% the mean color defined by 'color' (see above)
%
% AFQ_RenderFibers(fg,'subdivs',Nsubdivs) - The number of faces to create
% for each fiber node.  For example Nsubdivs = 6 would mean that each
% section of the fiber is hexagonal.  Increase Nsubdivs to make the fibers
% look like tubes. Default is Nsubdivs = 25
%
% AFQ_RenderFibers(fg,'tubes', [1]) - Whether fibers should be rendered as
% tubes [1] or lines [0]. Lines are faster tubes are prettier. Tubes are
% default
%
% AFQ_RenderFibers(fg,'newfig', [1]) - Check if rendering should be in new
% window [1] or added to an old window [0]. Default is new figure.
%
% AFQ_RenderFibers(fg, 'subplot', [row col num]) - Render the fiber in a
% subplot. This is like using the command subplot(row,col,num) to designate
% the position of the figure in the subplot.
%
% AFQ_RenderFibers(fg,'numfibers', numfibers) - Only render as many fibers
% as are defined in numfibers. This many fibers will randomly be chosen
% from the fiber group.  This is useful to save time when rendering large
% groups.
%
% AFQ_RenderFibers(fg,'alpha',alpha) - Set the transparency of the fibers.
% Alpha should be a value between 0 (transparent) and 1 (opaque).
%
% Outputs
% lightH    - Handle to the lighting object. This allows you to change the
%             lighting (e.g. camlight(lightH,'left');)
% fiberMesh - The fibers represented as a surface (see surf.m)
% fvc       - Face, vertex, color representation of fibers. This is
%             compatible with matlab's patch function and is a typical mesh
%             format used by other software.
%
% Example:
%
% fg = dtiReadFibers('Left_CST.mat');
% dt = dtiLoadDt6('dt6.mat');
% AFQ_RenderFibers(fg, 'dt', dt, 'radius', [.7 5], 'jittercolor', .1);
%
% Written by Jason D Yeatman, April 2012.
%
% Copyright Vista Team 2012

%% Argument checking

% Check to make sure the fiber group isn't empty
if isempty(fg.fibers) || length(fg.fibers) == 0
    fprintf('Fiber group is empty: %s\n',fg.name);
    return
end
% Check if a dt6 file was input
if sum(strcmpi('dt',varargin)) > 0
    dt = varargin{find(strcmpi('dt',varargin))+1};
else
    dt = [];
end

% Check if rois were input.  There should be 2 rois after the argument rois
if sum(strcmpi('rois',varargin)) > 0
    roi1 = varargin{find(strcmpi('rois',varargin))+1};
    roi2 = varargin{find(strcmpi('rois',varargin))+2};
    % Check if the rois are the propper format.  If not, do not use them
    if ~isstruct(roi1)
        warning('roi1 is not the propper format')
        roi1 = [];
        roi2 = [];
    elseif ~isstruct(roi2)
        warning('roi2 is not the propper format')
        roi1 = [];
        roi2 = [];
    end
else
    roi1 = [];
    roi2 = [];
end

% check if the camera angle is defined
if sum(strcmpi('camera',varargin)) > 0
    camera = varargin{find(strcmpi('camera',varargin))+1};
    if ischar(camera)
        switch(camera)
            case {'sagittal' 'leftsagittal' 'leftsag'}
                camera = [270 0];
                lightPosition = [-60 0 0];
            case {'rightsagittal' 'rightsag'}
                camera = [90 0];
                lightPosition = [60 0 0];
            case 'axial'
                camera = [0 90];
                lightPosition = [-10 10 80];
            case 'coronal'
                camera = [0 0];
                lightPosition = [-60 0 0];
            case 'uncinate'
                camera = [-85 40];
                lightPosition = [0 0 -100];
            case 'ifof'
                camera = [280 10];
                lightPosition = [-60 -10 -20];
            case 'callosum'
                camera = [180 90];
                lightPosition = [60 0 0];
        end
    else
        % default light position
        lightPosition = [-60 0 0];
    end
else
    camera = [270 0]; %default camera angle is looking at the sagital plane
    lightPosition = [-60 0 0];
end

% if color is not defined use gray
if sum(strcmpi('color',varargin)) > 0
    color = varargin{find(strcmpi('color',varargin))+1};
else
    color = [0.7 0.7 0.7]; %default color for fiber group is gray
end

% To make each fiber stand out more we can jitter the darkness of each
% fiber. The amount of change in darkness is randomly sampled from a
% uniform distributions with range = -js to js.
if sum(strcmpi('jittershading',varargin)) > 0
    js = varargin{find(strcmpi('jittershading',varargin))+1};
else
    js = .2; %default amount of jitter in the coloring
end

% To make each fiber stand out more we can jitter the color of each fiber.
% The amount of change in color is randomly sampled from a uniform
% distributions with range = -jf to jf.
if sum(strcmpi('jittercolor',varargin)) > 0
    jf = varargin{find(strcmpi('jittercolor',varargin))+1};
else
    jf = 0; %default amount of jitter in the coloring
end

% Number of subdivisions per fiber node
if sum(strcmpi('subdivs',varargin)) > 0
    subdivs = varargin{find(strcmpi('subdivs',varargin))+1};
else
    subdivs = 20;
end

% color range for tract profile
if sum(strcmpi('crange',varargin)) > 0
    crange = varargin{find(strcmpi('crange',varargin))+1};
else
    crange = [.3 .6]; %default color range
end

% color map for tract profile
if sum(strcmpi('cmap',varargin)) > 0
    cmap = varargin{find(strcmpi('cmap',varargin))+1};
else
    cmap = 'jet'; %default color map
end

% If a tract proflie is provided use it.  Otherwise if a dt6 file is input
% then compute the tract profile.  If neither is provided then do nothing.
if sum(strcmpi('tractprofile',varargin)) > 0
    TP = varargin{find(strcmpi('tractprofile',varargin))+1};
    computeTP = false;
elseif exist('dt','var') && ~isempty(dt)
    computeTP = true;
    TP = [];
else
    TP = [];
    computeTP = false;
end

% Check what value should be plotted on the tract profile
if sum(strcmpi('val',varargin)) > 0
    valname = varargin{find(strcmpi('val',varargin))+1};
else
    valname = 'fa';
end
% Check if the user has defined a radius for the fiber tubes and for the
% tract profile.  Within the radius variable the first number is for the
% fibers and the second number is for the tract profile.
if sum(strcmpi('radius',varargin)) > 0
    rFib = varargin{find(strcmpi('radius',varargin))+1};
    if length(rFib)==2
        rTP = rFib(2);
        rFib   = rFib(1);
    else
        % Set the radius of the tract profile to the default if the user
        % only set the radius of the fibers
        rTP = 5;
    end
else
    rFib = 0.7;
    rTP = 5;
end

% Check if fibers should be rendered as tubes or plotted as lines
if sum(strcmpi('tubes',varargin)) > 0
    tubes = varargin{find(strcmpi('tubes',varargin))+1};
else
    % default is tubes
    tubes = true;
end

% Check if rendering should be in new window or added to an old window
if sum(strcmpi('newfig',varargin)) > 0
    newfig = varargin{find(strcmpi('newfig',varargin))+1};
else
    % default is to render in a new figure window
    newfig = 1;
end

% Check if the figure should be in a subplot
if sum(strcmpi('subplot',varargin)) > 0
    splotnum = varargin{find(strcmpi('subplot',varargin))+1};
end

% Randomly choose the desired number of fibers to be rendered.
if sum(strcmpi('numfibers',varargin)) > 0
    numfib = varargin{find(strcmpi('numfibers',varargin))+1};
    % check if there are more fibers than the number specified
    if length(fg.fibers) > numfib
        % generate a random index of fibers
        fibindx = ceil(rand(numfib,1).*length(fg.fibers));
        % fibindx = randsample(length(fg.fibers),numfib);
        % retain only these fibers for the rendering
        fg.fibers = fg.fibers(fibindx);
        % if there fiber specific coloring was defined make sure to retain
        % the correct color values
        if iscell(color)
            color = color(fibindx);
        elseif size(color,1) > 1
            color = color(fibindx,:);
        end
    end
end

% Set the transparency of the fibers
if sum(strcmpi('alpha',varargin)) > 0
    alpha = varargin{find(strcmpi('alpha',varargin))+1};
else
    alpha = 1;
end
%% Loop over fibers and render them in 3-D

% Only render in a new figure window if desired (default)
if newfig == 1 && ~exist('splotnum','var')
    figure; hold('on');
elseif exist('splotnum','var') && ~isempty('splotnum')
    subplot(splotnum(1),splotnum(2),splotnum(3)); hold('on');
end

% Render each fiber as a tube (slow) if the parameter tubes == 1
if tubes == 1
    for ii = 1:length(fg.fibers)
        
        % X, Y and Z coordinates for the fiber
        coords = fg.fibers{ii};
        
        % Build the mesh of the fiber. The coloring of the mesh will depend
        % on whether color was defined for the full fiber group, for each
        % fiber, or for each point on each fiber
        if size(color,1) == 1 && size(color,2) == 3
            % If 1 color was defined than use that color for each fiber
            [X Y Z C] = AFQ_TubeFromCoords(coords, rFib, color, subdivs);
        elseif size(color,1) == length(fg.fibers) && size(color,2) == 3
            % If color is a N x 3 array with a different color defined for
            % each fiber than color fg.fibers{ii} the color defined by
            % color(ii,:)
            fColor = color(ii,:);
            [X, Y, Z, C,] = AFQ_TubeFromCoords(coords, rFib, fColor, subdivs);
        elseif iscell(color) && length(color) == length(fg.fibers)
            % If color is a cell array where each cell contains a color for
            % each point on the fiber than color each node fg.fibers{ii}
            % the colors defined in color{ii}
            fColor = color{ii};
            [X, Y, Z, C] = AFQ_TubeFromCoords(coords, rFib, fColor, subdivs);
        end
        
        % Add some random jitter to the shading so some fibers are slightly
        % lighter or darker
        C = C+[rand(1).*2 - 1].*js;
        % to jitter based on a Gaussian distribution
        % C = C+randn(1).*js;
        
        % Make sure that after adding the jitter the colors don't exceed the
        % range of 0 to 1
        C(C > 1) = 1;
        C(C < 0) = 0;
        
        % Add some random jitter to the fiber color.
        % jitter red channel
        C(:,:,1) = C(:,:,1)+[rand(1).*2 - 1].*jf;
        % jitter green channel
        C(:,:,2) = C(:,:,2)+[rand(1).*2 - 1].*jf;
        % jitter blue channel
        C(:,:,3) = C(:,:,3)+[rand(1).*2 - 1].*jf;
        % Make sure that after adding the jitter the colors don't exceed the
        % range of 0 to 1
        C(C > 1) = 1;
        C(C < 0) = 0;
        
        % Render fiber tubes and get the handle of this plot object
        h(ii) = surf(X,Y,Z,C,'facealpha',alpha);
        
        % Collect fiber coordinates to be returned
        if nargout > 1
            fiberMesh.X{ii} = X;
            fiberMesh.Y{ii} = Y;
            fiberMesh.Z{ii} = Z;
            fiberMesh.C{ii} = C;
        end
    end
    
    % Plot the fibers as lines (much faster than tubes) if the input tubes == 0
else
    % If only one color is provided then repeat it for each fiber
    if size(color,1) == 1
        color = repmat(color,length(fg.fibers),1);
    end
    for ii = 1:length(fg.fibers)
        % X, Y and Z coordinates for the fiber
        x = fg.fibers{ii}(1,:)';
        y = fg.fibers{ii}(2,:)';
        z = fg.fibers{ii}(3,:)';
        % Jitter color and shading
        C = color(ii,:) + [rand(1,3).*2 - 1].*jf+[rand(1).*2 - 1].*js;
        % Make sure that after adding the jitter the colors don't exceed the
        % range of 0 to 1
        C(C > 1) = 1;
        C(C < 0) = 0;
        % plot the fibers as lines
        plot3(x,y,z,'-','color',C);
        
        % Collect fiber coordinates to be returned
        if nargout > 1
            fiberMesh.X{ii} = x;
            fiberMesh.Y{ii} = y;
            fiberMesh.Z{ii} = z;
            fiberMesh.C{ii} = C;
        end
    end
end

%% Render the tract profile
if computeTP
    % compute tract profile
    [fa,md,rd,ad,cl,fgCore]=dtiComputeDiffusionPropertiesAlongFG(fg, dt, roi1, roi2, 100);
    % Create Tract Profile structure
    TP = AFQ_CreateTractProfile;
    % Set the desired values to the structure
    TP = AFQ_TractProfileSet(TP,'vals',valname,eval(valname));
    % Set the tract profile coordinates
    TP = AFQ_TractProfileSet(TP,'coordsacpc',fgCore.fibers{1});
end
if exist('TP','var') && ~isempty(TP)
    % Get the coordinates
    coords = AFQ_TractProfileGet(TP,'coordsacpc');
    % Get the values
    vals = AFQ_TractProfileGet(TP,'vals',valname);
    % Render the tract Profile
    AFQ_RenderTractProfile(coords,rTP,vals,30,cmap,crange);
end

% Render the ROIs if they were passed in
if exist('roi1','var')&&exist('roi2','var')&&~isempty(roi1)&&~isempty(roi2)
    AFQ_RenderRoi(roi1,[1 0 0]);
    AFQ_RenderRoi(roi2,[1 0 0]);
end
%% Set the lighting, axes, etc.

shading('interp');
% Only set window properties if the fibers were rendered in a new figure
% window
if newfig ==1
    coordsAll = horzcat(fg.fibers{:});
    mx = minmax(coordsAll(1,:));
    my = minmax(coordsAll(2,:));
    mz = minmax(coordsAll(3,:));
    axis([mx(1)-3 mx(2)+3 my(1)-3 my(2)+3 mz(1)-3 mz(2)+3],'equal');
    xlabel('X mm','fontname','times','fontsize',14);
    ylabel('Y mm','fontname','times','fontsize',14);
    zlabel('Z mm','fontname','times','fontsize',14);
    set(gca,'fontname','times','fontsize',14);
    set(gcf,'color',[1 1 1]);
    grid('on');
    if size(camera,2) == 2
        view(camera(1), camera(2));
    elseif size(camera,2) == 3
        campos(camera);
    end
    lightH = camlight('right');
    %set(lightH, 'position',lightPosition);
    lighting('gouraud');
    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit');
    fprintf('\nmesh can be rotated with arrow keys\n')
else
    lightH = [];
end

%% Return the mesh of the fibrs as a patch object
% if tubes == 1
%     for ii = 1:length(h)
%         fvc(ii) = surf2patch(h(ii));
%     end
% end
