function lightH = AFQ_RenderFibers(fg,varargin)
% Render a fiber group and tract profile in 3-D
%
%   AFQ_RenderFibers(fg,'PropertyName',PropertyValue, ...)
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
% AFQ_RenderFibers(fg,'tractprofile',TP) - A tract profile can be passed in
% if it has been precomputed for the fiber tract. TP is a structure with
% fields TP.vals and TP.coords that define the values to be plotted and the
% coordinates of the tract profile.  See AFQ_CreateTractProfile.
%
% AFQ_RenderFibers(fg,'camera', view) - Render the fiber group and view
% from a specific plane.  The options for view are 'sagittal', 'coronal' or
% 'axial'. The default is sagittal. View can also be defined with a 1 by 2
% vector denoting the azimuth (horizontal rotation) and elevation (vertical
% rotation) in degrees of the camera with respect to the fibers. See the
% matlab view function
%
% AFQ_RenderFibers(fg,'color', rgbValues) - Render the fiber group in a
% specific color.  Color is a 1 by 3 vector of rgb values. The default is
% gray [0.7 0.7 0.7]. To do cyan for example rgbValues = [0 1 1].
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
            case 'sagittal'
                camera = [270 0];
                lightPosition = [-60 0 0];
            case 'axial'
                camera = [0 90];
                lightPosition = [-60 0 0];
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
%     if length(camera) ~= 2
%         camera = [270 0];
%         lightPosition = [-60 0 0];
%         warning('Camera angle needs 2 numbers. Set to default')
%     end
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
    subdivs = 25;
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
    rFib = 1;
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
%% Loop over fibers and render them in 3-D

% Only render in a new figure window if desired (default)
if newfig == 1
    figure; hold('on');
end

% Render each fiber as a tube (slow) if the parameter tubes == 1
if tubes == 1
    for ii = 1:length(fg.fibers)
        % X, Y and Z coordinates for the fiber
        x = fg.fibers{ii}(1,:)';
        y = fg.fibers{ii}(2,:)';
        z = fg.fibers{ii}(3,:)';
        
        % Initialize the variables for the mesh
        N = length(x);
        X=zeros(N,subdivs);
        Y=zeros(N,subdivs);
        Z=zeros(N,subdivs);
        theta=0:(2*pi/(subdivs-1)):(2*pi);
        
        % frame seems to work better than frenet
        [t,n,b]=frame(x,y,z,randn(1,3));
        %[t,n,b]=frenet(x,y,z);
        
        % set the radius of the tubes
        r=rFib*ones(N,1);
        
        % Build a mesh for the fibers
        for i=1:N
            X(i,:)=x(i) + r(i)*(n(i,1)*cos(theta) + b(i,1)*sin(theta));
            Y(i,:)=y(i) + r(i)*(n(i,2)*cos(theta) + b(i,2)*sin(theta));
            Z(i,:)=z(i) + r(i)*(n(i,3)*cos(theta) + b(i,3)*sin(theta));
        end
        
        % Set the color for each point on the mesh
        C = ones(size(Z,1),size(Z,2),3);
        C(:,:,1) = C(:,:,1).*color(1);
        C(:,:,2) = C(:,:,2).*color(2);
        C(:,:,3) = C(:,:,3).*color(3);
        
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
        
        % Render fiber tubes
        surf(X,Y,Z,C);
    end
    
    % Plot the fibers as lines (much faster than tubes) if the input tubes == 0
else
    for ii = 1:length(fg.fibers)
        % X, Y and Z coordinates for the fiber
        x = fg.fibers{ii}(1,:)';
        y = fg.fibers{ii}(2,:)';
        z = fg.fibers{ii}(3,:)';
        % Jitter color and shading
        C = color + [rand(1,3).*2 - 1].*jf+[rand(1).*2 - 1].*js;
        % Make sure that after adding the jitter the colors don't exceed the
        % range of 0 to 1
        C(C > 1) = 1;
        C(C < 0) = 0;
        % plot the fibers as lines
        plot3(x,y,z,'-','color',C);
    end
end

%% Render the tract profile
if computeTP
    % compute tract profile
    [fa,md,rd,ad,cl,fgCore]=dtiComputeDiffusionPropertiesAlongFG(fg, dt, roi1, roi2, 100);
    TP = AFQ_CreateTractProfile;
    TP.vals.fa = fa;
    TP.coords  = fgCore.fibers{1};
end
if exist('TP','var') && ~isempty(TP)
    AFQ_RenderTractProfile(TP.coords,rTP,TP.vals.fa,30,cmap,crange);
end
%% Set the lighting, axes, etc.
coordsAll = horzcat(fg.fibers{:});
mx = minmax(coordsAll(1,:));
my = minmax(coordsAll(2,:));
mz = minmax(coordsAll(3,:));
axis([mx(1)-3 mx(2)+3 my(1)-3 my(2)+3 mz(1)-3 mz(2)+3],'equal');
% Only set window properties if the fibers were rendered in a new figure
% window
if newfig ==1
    xlabel('X mm','fontname','times','fontsize',14);
    ylabel('Y mm','fontname','times','fontsize',14);
    zlabel('Z mm','fontname','times','fontsize',14);
    set(gca,'fontname','times','fontsize',14);
    set(gcf,'color',[1 1 1]);
    shading('interp');
    grid('on');
    if size(camera,2) == 2
        view(camera(1), camera(2));
    elseif size(camera,2) == 3
        campos(camera)
    end
    lightH = camlight;
    set(lightH, 'position',lightPosition);
    lighting('gouraud');
    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit');
    fprintf('\nmesh can be rotated with arrow keys\n')
end
