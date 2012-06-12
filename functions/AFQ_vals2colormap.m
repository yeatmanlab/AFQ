function rgb = AFQ_vals2colormap(vals, colormap, crange)
% Take in a vector of N values and return and return a Nx3 matrix of RGB
% values associated with a given colormap
%
% rgb = AFQ_vals2colormap(vals, [colormap = 'jet'], [crange])
%
% Inputs:
% vals     = A vector of values to map to a colormap
% colormap = A matlab colormap. Examples: colormap = 'autumn'; 
%            colormap = 'jet'; colormap = 'hot';
% crange   = The values to map to the minimum and maximum of the colormap.
%            Defualts to the full range of values in vals.
%
% Outputs:
% rgb      = Nx3 matrix of rgb values mapping each value in vals to the
%            corresponding colormap value
%
% Example:
% vals = rand(1,100);
% rgb = AFQ_vals2colormap(vals, 'hot');
%
% Copyright Jason D. Yeatman, June 2012

if ~exist('colormap','var') || isempty(colormap)
    colormap = 'jet';
end
if ~exist('crange','var') || isempty(crange)
    crange = [min(vals) max(vals)];
end
% Generate the colormap
cmap = eval([colormap '(256)']);
% Normalize the values to be between 1 and 256
vals(vals < crange(1)) = crange(1);
vals(vals > crange(2)) = crange(2);
valsN = round(((vals - crange(1)) ./ diff(crange)) .* 255)+1;
% Convert any nans to ones
valsN(isnan(valsN)) = 1;
% Convert the normalized values to the RGB values of the colormap
rgb = cmap(valsN, :);

return
