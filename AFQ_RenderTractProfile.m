function [lightH] = AFQ_RenderTractProfile(coords, radius, color, subdivs, cmap, crange)

% AFQ_RenderTractProfile - plots a tube r along the space curve x,y,z.
%
% AFQ_RenderTractProfile(x,y,z) plots the basic tube with radius 1
% AFQ_RenderTractProfile(x,y,z,r) plots the basic tube with variable radius r (either
% a vector or a value) AFQ_RenderTractProfile(x,y,z,r,v) plots the basic tube with
% coloring dependent on the values in the vector v
% AFQ_RenderTractProfile(x,y,z,r,v,s) plots the tube with s tangential subdivisions
% (default is 6)
%
% [X,Y,Z]=AFQ_RenderTractProfile(x,y,z) returns [Nx3] matrices suitable for mesh or
% surf
%
% Note that the tube may pinch at points where the normal and binormal
% misbehaves. It is suitable for general space curves, not ones that
% contain straight sections. Normally the tube is calculated using the
% Frenet frame, making the tube minimally twisted except at inflexion
% points.
%
% To deal with this problem there is an alternative frame:
% AFQ_RenderTractProfile(x,y,z,r,v,s,vec) calculates the tube by setting the normal
% to the cross product of the tangent and the vector vec. If it is chosen
% so that it is always far from the tangent vector the frame will not twist
% unduly
%
% Example:
%
%  t=0:(2*pi/100):(2*pi); x=cos(t*2).*(2+sin(t*3)*.3);
%  y=sin(t*2).*(2+sin(t*3)*.3); z=cos(t*3)*.3;
%  AFQ_RenderTractProfile(x,y,z,0.14*sin(t*5)+.29,t,10)
%
% Derived from code by Anders Sandber, asa@nada.kth.se, 2005
% Modified by Jason D Yeatman, April 2012. 
%
% Copyright Vista Team 2012

%% Argument checking

% if the colormap is not defined set to the defaul, jet.
if ~exist('cmap','var') || isempty(cmap)
    cmap = 'jet'
end
% if the color range is not defined set to the defaul, jet.
if ~exist('crange','var') || isempty(crange)
    crange = [.3 .6]; % defualt color range
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
lightH = camlight('headlight');
axis('equal');
% set the axis and the color range
axis([min(x)-r(1) max(x)+r(1) min(y)-r(1) max(y)+r(1)...
    min(z)-r(1) max(z)+r(1) crange(1) crange(2)]);
% add a colorbar
colorbar;

% set(gca,'cameraposition',[min(x)-r(1) 0 0],'cameratarget',[0 0 0]);
% xlabel('X mm'); ylabel('Y mm'); zlabel('Z mm')
% camlight('headlight');

return
