function [lightH] = AFQ_RenderTractProfile(coords, radius, color, subdivs, cmap, crange, newfig)

% Render Tract Profile as a tube colored based on the Tract Profile values
%
%   [lightH] = AFQ_RenderTractProfile(coords, radius, color, subdivs, cmap, crange, [newfig=0])
%
% Inputs:
% 
% coords  = N by 3 matrix of x, y, z coordinates of the tract profile
% radius  = Desired radius of the rendered tube
% color   = A vector of N values that will denote the color at each point on
%           the Tract Profile.
% subdivs = Number of faces for each section of the tube.  More means the
%           tube is more circular. Less subdivs makes tube blocky
% cmap    = Desired colormap. eg. 'jet' 'hot' 'autumn' 'winter' 
% crange  = Maximum and minimum values of the colormap
% newfig  = whether or not to open a new figure window
%
% Outputs:
%
% lightH  = A handle for the lighting so it can be adjusted later
%
% Example:
%
% core    = dtiComputeSuperFiberRepresentation(fg);
% coords  = core.fibers{1};
% radisu  = 5;
% color   = tiComputeDiffusionPropertiesAlongFG(fg);
% subdivs = 20;
% cmap    = 'jet';
% crange  = [0.3 0.6];
% [lightH] = AFQ_RenderTractProfile(coords, radius, color, subdivs, cmap, crange)
% 
% Derived from code by Anders Sandber, asa@nada.kth.se, 2005
% Modified by Jason D Yeatman, April 2012. 
%
% Copyright Vista Team 2012

%% Argument checking

% if the colormap is not defined set to the defaul, jet.
if ~exist('cmap','var') || isempty(cmap)
    cmap = 'jet';
end
% if the color range is not defined set to the defaul, jet.
if ~exist('crange','var') || isempty(crange)
    crange = [.3 .6]; % defualt color range
end
% default is no new figure window
if ~exist('newfig','var') || isempty(newfig)
    newfig = 0;
elseif newfig == 1
    figure;
end
% Compute number of coordinates
N=size(coords,1);
% transpose if necesary
if (N == 3)
    coords = coords';
    N = size(coords,1);
end
% divide into x, y, z
x = coords(:,1);
y = coords(:,2);
z = coords(:,3);
% Set radius of tract profile
if ~exist('radius','var') || isempty(radius)
    r=x*0+1;
else
    r = radius;
    if (size(r,1)==1 & size(r,2)==1)
        r=r*ones(N,1);
    end
end
% Set the number of subdivisions for the tube
if ~exist('subdivs','var') || isempty(subdivs)
    subdivs=25;
elseif length(subdivs) > 1
    error('subdivs must be a scaler')
end
% frame seems to work better than frenet.  With frenet sometimes the tube
% will pinch at points.  A random vector seems to work well.
vec = randn(1,3);
[t,n,b]=frame(x,y,z,vec);
%[t,n,b]=frenet(x,y,z);


%% Build the mesh
X=zeros(N,subdivs);
Y=zeros(N,subdivs);
Z=zeros(N,subdivs);
theta=0:(2*pi/(subdivs-1)):(2*pi);
for i=1:N
    X(i,:)=x(i) + r(i)*(n(i,1)*cos(theta) + b(i,1)*sin(theta));
    Y(i,:)=y(i) + r(i)*(n(i,2)*cos(theta) + b(i,2)*sin(theta));
    Z(i,:)=z(i) + r(i)*(n(i,3)*cos(theta) + b(i,3)*sin(theta));
end

%% Render the mesh
if exist('color','var') && ~isempty(color)
    % color based on the specified color vector
    V = color;
    if (size(V,1)==1)
        V=V';
    end
    V = V*ones(1,subdivs);
    surf(X,Y,Z,V);
else
    % if no color vector is suplied then render with the defaults
    surf(X,Y,Z);
end

%% Setting the lighting etc.
colormap(cmap) 
shading('interp');
lighting('gouraud');
% set the axis and the color range
% axis([min(x)-r(1) max(x)+r(1) min(y)-r(1) max(y)+r(1)...
%     min(z)-r(1) max(z)+r(1) crange(1) crange(2)]);
caxis(crange);
% add a colorbar
colorbar;

% Set figure window properties if it is a new window
if newfig == 1
    lightH = camlight('headlight');
    axis('image');
end

% set(gca,'cameraposition',[min(x)-r(1) 0 0],'cameratarget',[0 0 0]);
% xlabel('X mm'); ylabel('Y mm'); zlabel('Z mm')
% camlight('headlight');

return
