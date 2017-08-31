function rgb = vals2colormap(vals, colormap, crange)
% Take in a vector of N values and return and return a Nx3 matrix of RGB
% values associated with a given colormap
%
% rgb = AFQ_vals2colormap(vals, [colormap = 'jet'], [crange])
%
% Inputs:
% vals     = A vector of values to map to a colormap or a cell array of
%            vectors of values
% colormap = A matlab colormap. Examples: colormap = 'autumn';
%            colormap = 'jet'; colormap = 'hot';
% crange   = The values to map to the minimum and maximum of the colormap.
%            Defualts to the full range of values in vals.
%
% Outputs:
% rgb      = Nx3 matrix of rgb values mapping each value in vals to the
%            corresponding rgb colors.  If vals is a cell array then rgb
%            will be a cell array of the same length
%
% Example:
% vals = rand(1,100);
% rgb = AFQ_vals2colormap(vals, 'hot');
%
% Copyright Jason D. Yeatman, June 2012

if ~exist('colormap','var') || isempty(colormap)
    colormap = 'jet';
end
% Generate the colormap
if ischar(colormap)
    cmap = eval([colormap '(256)']);
else
    % If the colors were provided calculate the number of unique colors
    ncolors = size(colormap,1);
    rp = floor(256./ncolors);
    rm = 256 - rp.*ncolors;
    cmap = [];
    % Repeat each color an equal number of times
    for ii = 1:ncolors
       cmap = vertcat(cmap,repmat(colormap(ii,:),rp,1));
    end
    cmap = vertcat(repmat(colormap(1,:),round(rm./2),1),cmap,repmat(colormap(end,:),round(rm./2),1));
end
%
if ~iscell(vals)
    if ~exist('crange','var') || isempty(crange)
        crange = [min(vals) max(vals)];
    end
    % Normalize the values to be between 1 and 256
    vals(vals < crange(1)) = crange(1);
    vals(vals > crange(2)) = crange(2);
    valsN = round(((vals - crange(1)) ./ diff(crange)) .* 255)+1;
    % Convert any nans to ones
    valsN(isnan(valsN)) = 1;
    % Convert the normalized values to the RGB values of the colormap
    rgb = cmap(valsN, :);
elseif iscell(vals)
    if ~exist('crange','var') || isempty(crange)
        if size(vals{1},1) > size(vals{1},2)
            crange = [min(vertcat(vals{:})) max(vertcat(vals{:}))];
        else
            crange = [min(horzcat(vals{:})) max(horzcat(vals{:}))];
        end
    end
    for ii = 1:length(vals)
        % Normalize the values to be between 1 and 256 for cell ii
        valsN = vals{ii};
        valsN(valsN < crange(1)) = crange(1);
        valsN(valsN > crange(2)) = crange(2);
        valsN = round(((valsN - crange(1)) ./ diff(crange)) .* 255)+1;
        % Convert any nans to ones
        valsN(isnan(valsN)) = 1;
        % Convert the normalized values to the RGB values of the colormap
        rgb{ii} = cmap(valsN, :);
    end
end
return
