function h = AFQ_AddImageTo3dPlot(nifti, slice, cmap, rescale, alpha, imgClipRange, varargin)
% Add a slice from a nifti image to a 3D plot
%
% h = AFQ_AddImageTo3dPlot(nifti, slice, [cmap], [rescale], alpha, imgClipRange, varargin)
%
% Inputs:
%
% nifti   = A 3D nifti image. It is ok if the image has translation
%           information in the header, however there may be problems if
%           there are rotations saved into the header. Typically we do not
%           save rotations into the header so everything should be fine.
%
% slice   = Designate which slice in acpc (millimeter) coordinates should
%           be added to the plot. slice should be a 1x3 vector where each
%           coordinate is either a 0 or nan except the coordinate of the
%           slice to be rendered. For example to render the coronal slice 5
%           mm anterior to the anterior commissure: slice = [5 0 0]. For an
%           axial slice 10mm below the acpc plane: slice = [0 0 -10]. For a
%           sagittal slice 25mm to the left of the mid-sagittal plane:
%           slice = [-25 0 0]
%
% cmap    = Define a colormap for the values in the image.  The defualt is
%           cmap = 'gray'; Other options: cmap = 'jet'; cmap = 'hot'; etc.
%
% rescale = If the image is being added to an existing axis the default is
%           to rescale the axis so the full image can be seen 
%           [rescale = 1]. To preserve the current axis properties: 
%           rescale = 0
% 
% alpha   = Alpha of the image. Value from 0 to 1 that sets how 
%           transparent the image is.
%
% imgClipRange = The image will be scaled to be between these two values.
%                If this is left blank then we will try to guess a good
%                range of values to make it look nice.
%
% Aditional options:
% 
% To add an overlay heatmap to the rendering. For example BOLD activation.
% Add the aditional argument 'overlay', followed by a nifti image and
% threshold.
% AFQ_AddImageTo3dPlot(nifti, slice, [cmap], [rescale], alpha,'overlay',nifti,threshold)
% 
% Output:
%
% h       = Handle for the object in the figure.  You can remove the slice
%           from the figure window with: delete(h);
%
% Example:
%
% % Load a fiber group
% fg = dtiReadFibers('LeftArcuate.mat');
% % Render the fibers as red tubes.
% AFQ_RenderFibers(fg, 'color', [1 0 0]);
% % Load a nifti image
% nifti = readFileNifti('t1.nii.gz');
% % Add a sagittal slice to the figure;
% h = AFQ_AddImageTo3dPlot(nifti, [-20 0 0])
%
% Copyright Jason D Yeatman June 2012

%% Check arguments
if ~exist('cmap','var') || isempty(cmap)
    cmap = gray(256);
elseif ischar(cmap)
    cmap = eval([cmap '(256)']);
end
if ~exist('rescale','var') || isempty(rescale)
    rescale = 1;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end
if ~exist('imgClipRange','var') || isempty(imgClipRange)
    imgClipRange = [];
end
% Make sure only one plane was designated for plotting
plane = ~isnan(slice);
if sum(plane) > 1
    plane = slice ~= 0;
end
if sum(plane) > 1
    error('Only one slice can be plotted. Other values should be nan or 0')
end
% Replace nans with zeros
slice = slice .* plane;
% Make sure slice is a row vector
if size(slice,2) ~= 3
    slice = slice';
end
%% Create a color image from the slice in the 3d volume
[colorimg min_x max_x min_y max_y min_z max_z] = ...
    MakeImageFromVolume(nifti,slice,plane,cmap,imgClipRange);

%% Create overlay image
% Check if an overlay image was passed in
ov = find(strcmpi('overlay',varargin));
if sum(ov) == 1
    % Get the overlay image volume
    overlayIm = varargin{ov+1};
    % Get the overlay threshold
    overlayThresh = varargin{ov+2};
    % Get the colormap if one was defined
    if length(varargin) > ov+2
        overlayCmap = varargin{ov+3}
    else
        overlayCmap = autumn(256);
    end
    % Create a color image
    OverColorImg = MakeImageFromVolume(overlayIm,slice,plane,overlayCmap,overlayThresh);
    % Combine the overlay and the background images
    nonan = find(~isnan(OverColorImg));
    colorimg(nonan) = OverColorImg(nonan);
end

%% Plot the image
% The call to surface will be slightly different depending on whether it is
% an axial, sagittal or coronal image.
if find(plane) == 1
    % Add the image to the current 3D plot
    z_dims = [min_z max_z; min_z max_z];
    h = surf([min_x max_x],[min_y max_y],z_dims,...
        colorimg,'facecolor','texture','ambientstrength',1,'facealpha',alpha);
elseif find(plane) == 2
    % The image dimensions must be permuted to match what is expected in
    % surf
    colorimg=permute(colorimg,[2 1 3]);
    % Add the image to the current 3D plot
    z_dims = [min_z min_z;   max_z max_z];
    h = surf([min_x max_x],[min_y max_y],z_dims,...
        colorimg,'facecolor','texture','ambientstrength',1,'facealpha',alpha);
elseif find(plane) == 3
    % The image dimensions must be permuted to match what is expected in
    % surf
    colorimg=permute(colorimg,[2 1 3]);
    % Add the image to the current 3D plot
    h =surf([min_x max_x],[min_y max_y],repmat(min_z, [2 2]),...
        colorimg,'facecolor','texture','ambientstrength',1,'facealpha',alpha);
end

% Rescale the axis to fit the image
if rescale == 1
    % Old axis values
    ax = axis;
    % New axis scaling
    ax2 = [min_x max_x min_y max_y min_z max_z];
    % Only change the dimensions that are needed to see the full image
    if find(plane) == 1
        ax2(1:2) = ax(1:2);
    elseif find(plane) == 2
        ax2(3:4) = ax(3:4);
    elseif find(plane) == 3;
            a = get(gca);
            if isfield(a,'ZLim')
                ax2(5:6) = a.ZLim;
            end
    end
    % Rescale axis
    axis(ax2);
end

% Turn hold on in case other features are added to the rendering
hold on;

return

function [colorimg min_x max_x min_y max_y min_z max_z] =...
     MakeImageFromVolume(nifti,slice,plane,cmap,thresh)
% Get a slice from a nifti volume, resample it to 1mm isotropic resolution
% and make a rgb image from it

% First Format the image with propper transforms etc.
% The variable slice is given in acpc coordinates.  Transform acpc to image
% coordinates
imgXform = nifti.qto_ijk;
ImCoords = round(imgXform * [slice 1]');
imIndx = ImCoords(find(slice));

% Pull the desired slice out of the 3d image volume
if find(plane) == 1
    image = squeeze(nifti.data(imIndx,:,:));
    % Define the minimum and maximum coordinates of the image in each
    % dimension in acpc mm space.
    min_x = slice(1); max_x = min_x;
    [y z] = size(image);
    max_corner = inv(imgXform) * [imIndx y z 1]';
    max_y = max_corner(2); max_z = max_corner(3);
    min_corner = inv(imgXform) * [imIndx 0 0 1]';
    min_y = min_corner(2); min_z = min_corner(3);
elseif find(plane) == 2
    image = squeeze(nifti.data(:,imIndx,:));
    % Define the minimum and maximum coordinates of the image in each
    % dimension in acpc mm space.
    min_y = slice(2); max_y = min_y;
    [x z] = size(image);
    max_corner = inv(imgXform) * [x imIndx z 1]';
    max_x = max_corner(1); max_z = max_corner(3);
    min_corner = inv(imgXform) * [0 imIndx 0 1]';
    min_x = min_corner(1); min_z = min_corner(3);
else
    image = squeeze(nifti.data(:,:,imIndx));
    % Define the minimum and maximum coordinates of the image in each
    % dimension in acpc mm space.
    min_z = slice(3); max_z = min_z;
    [x y] = size(image);
    max_corner = inv(imgXform) * [x y imIndx 1]';
    max_x = max_corner(1); max_y = max_corner(2);
    min_corner = inv(imgXform) * [0 0 imIndx 1]';
    min_x = min_corner(1); min_y = min_corner(2);
end

% Resample the image to 1mm resolution. The necessary scale factor is
% stored in the image xform.
scale = diag(nifti.qto_xyz); scale = scale(1:3);

% The new dimensions will be the old dimensions multiplied by the scale
% factors for the plane
oldDim = size(image);
newDim = oldDim .* scale(find(plane == 0))';
% Resize the image
image = double(imresize(image,newDim));

if ~exist('thresh','var') || isempty(thresh)
    % Scale and clip the image values so that the lowest 1% of the values are
    % zeroed, the top 1% are maxed and the range is 0 to 255
    image = uint8(mrAnatHistogramClip(image,.01,.99,1).* 255);
else
    % If a threshold was defined then change values below that threshold to nan
    nanindx = image < thresh(1);
    image(nanindx) = thresh(1);
    if length(thresh) == 2
        % If a maximum threshold was defined then clip values above that
        image(image > thresh(2)) = thresh(2);
    end
    % Scale the image based on the threshold values
    imMin = min(image(:));
    image = image - imMin;
    imMax = max(image(:));
    image = uint8((image./imMax).*255);
end

% Convert the scaler image to an RGB image based on the chosen colormap
colorimg = ind2rgb(image,cmap);

% If a threshold was defined then change values below that threshold to nan
if exist('thresh','var') && ~isempty(thresh)
    nanindx3 = repmat(nanindx,[1 1 3]);
    colorimg(nanindx3) = nan;
end


