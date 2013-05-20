function params = CreateParamsStruct(c)
% Create a parameters structure from a cell array.
%
% params = CreateParamsStruct(c)
%
% Take a cell array and convert it to a structure where each field in the
% structure is named based on the odd numbered cells of the array and the
% values in each field come from the even numbered cells in the array. This
% function is useful for taking 'param', val, pairs from varargin and
% putting them in a structure

if ~isparams(c) && iscell(c)
    % If the cell array only contains one cell check if it is a params
    % structure and if so return it
    if length(c) == 1 && isparams(c{1})
        params = c{1};
    elseif length(c) == 1 && ~isparams(c{1})
        params = c;
    else
        params = struct;
        params.type = 'params version 1';
        for ii = 1:2:length(c)
            params.(c{ii}) = c{ii+1};
        end
    end
end