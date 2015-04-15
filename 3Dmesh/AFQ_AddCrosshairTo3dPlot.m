function h = AFQ_AddCrosshairTo3dPlot(coord, cross_size)
% Add crosshairs to a 3d plot
%
% h = AFQ_AddCrosshairTo3dPlot(coord, cross_size)
%
% Inputs:
%
% coord      = The location [x y z] for the crosshair
% cross_size = The length and radius of the arms of the crosshair.
%              [length radius]. For example [5 1] would mean each arm of
%              the crosshair is 5mm long with a radius of 1mm/
% Outputs:
%
% h          = Figure handles for the crosshair. h(1) is the x tube, h(2)
%              is the y tube and h(3) is the z tube.
%
% Copyright Jason D Yeatman, June 2012

% number of faces for each section of the tube
subdivs = 20;
% radius and length of tubes
L = cross_size(1);
R = cross_size(2);
% Coordinates for the tubes. First the X direction of the crosshair
x_x = [coord(1) - L : coord(1) + L];
x_y = repmat(coord(2), size(x_x,1), length(x_x));
x_z = repmat(coord(3), size(x_x,1), length(x_x));
% Y diretion of the crosshair
y_y = [coord(2) - L : coord(2) + L];
y_x = repmat(coord(1), length(x_x), length(x_x));
y_z = repmat(coord(3), length(x_x), length(x_x));
% Z direction of the crosshair
z_z = [coord(3) - L : coord(3) + L];
z_x = repmat(coord(1), length(x_x), length(x_x));
z_y = repmat(coord(2), length(x_x), length(x_x));
% Build and render meshes of each tube
h(1) = surf(AFQ_TubeFromCoords([x_x' x_y' x_z'],R,[1 0 0],subdivs));
h(2) = surf(AFQ_TubeFromCoords([y_x' y_y' y_z'],R,[0 1 0],subdivs));
h(3) = surf(AFQ_TubeFromCoords([z_x' z_y' z_z'],R,[0 0 1],subdivs));

% Turn hold on in case other features are added to the rendering
hold on;

return
