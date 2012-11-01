function [X Y Z C] = AFQ_TubeFromCoords(coords, r, color, subdivs)
% Build a mesh of a tube to be rendered with surf.m from a coordinate array
%    [X Y Z C] = AFQ_TubeFromCoords(coords, r, color, subdivs)
%
% Inputs:
%
% coords  = An array of x, y, z coordinates to be rendered as a tube
% r       = Radius of the tube
% color   = The color for the tube. Either a 1x3 vector containing the RGB
%           values for the surface of the tube or a Nx3 vector where N is
%           the number of coordinates.  In this case each coordinate will
%           have its own color.
% subdivs = The number of faces for each section of the tube
%
% Outputs:
%
% X       = X coordinates of the mesh
% Y       = Y coordinates of the mesh
% Z       = Z coordinates of the mesh
% C       = The color at each point on the mesh
%
% Example:
% % The coordinates of 1 fiber
% coords = fg.fibers{1};
% % Have the color change randomely along the fiber
% color = rand(size(coords,1),3);
% [X Y Z C] = AFQ_TubeFromCoords(coords, 1, color, 20)
% surf(X, Y, Z, C);
%
% Copyright Jason D Yeatman

if size(coords,1) == 3
    coords = coords';
end
% X, Y and Z coordinates for the tube
x = coords(:,1);
y = coords(:,2);
z = coords(:,3);
% Initialize the variables for the mesh
N = length(x);
X=zeros(N,subdivs);
Y=zeros(N,subdivs);
Z=zeros(N,subdivs);
theta=0:(2*pi/(subdivs-1)):(2*pi);

% frame seems to work better than frenet
[t,n,b]=frame(x,y,z,randn(1,3));
%[t,n,b]=frenet(x,y,z);

% set the radius of the tube
r=r*ones(N,1);

% Build a mesh for the fiber
for i=1:N
    X(i,:)=x(i) + r(i)*(n(i,1)*cos(theta) + b(i,1)*sin(theta));
    Y(i,:)=y(i) + r(i)*(n(i,2)*cos(theta) + b(i,2)*sin(theta));
    Z(i,:)=z(i) + r(i)*(n(i,3)*cos(theta) + b(i,3)*sin(theta));
end

% Set the color for each point on the mesh
if size(color,1) == 1 && size(color,2) == 3
    % If 1 color was defined than use that color for the fiber
    C = ones(size(Z,1),size(Z,2),3);
    C(:,:,1) = C(:,:,1).*color(1);
    C(:,:,2) = C(:,:,2).*color(2);
    C(:,:,3) = C(:,:,3).*color(3);
elseif size(color,1) == size(coords,1)
    % If color has a row for each fiber than each fiber node will have its
    % own color
    fcolor = reshape(color,[N 1 3]);
    C = repmat(fcolor,[1 subdivs 1]);
end