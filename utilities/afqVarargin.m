function yourStruct = afqVarargin(yourStruct,varArgIn)
%
%  fucntion yourStruct = afqVarargin(yourStruct,varArgIn)
%
%  This funciton is a wrapper/loop for varargin that takes as input a
%  structure and 'varargin' from the calling function and loops over
%  varargin to set the input structure's fields and check to make sure that
%  the field names in varargin are valid.
%
% INPUT:
%       yourStruct - a structure with fields
%       varArgIn   - varargin argument from the calling funciton.
%
% OUTPUT:
%       yourStruct - your original structure with fields set by varargin
%
%
% (C) Stanford University, VISTA 2011 [lmp]
%

%% Varargin

if ~isempty(varArgIn)
    for ii = 1:2:numel(varArgIn)-1
        if isfield(yourStruct,varArgIn{ii})
            % Check to make sure that the argument is formatted properly
            [param value] = afqCheckArgument(varArgIn{ii}, varArgIn{ii+1});
            yourStruct = setfield(yourStruct, param, value);
            clear param value
        else
            %Check for fields within fields
            fnames = fieldnames(yourStruct);
            for jj = 1:length(fnames)
                if isstruct(yourStruct.(fnames{jj})) && isfield(yourStruct.(fnames{jj}),varArgIn{ii})
                    % Check to make sure that the argument is formatted properly
                    [param value] = afqCheckArgument(varArgIn{ii}, varArgIn{ii+1});
                    yourStruct.(fnames{jj}) = setfield(yourStruct.(fnames{jj}), param, value);
                end
                %warning('"%s" is not a valid field name.\n',varArgIn{ii});
            end
        end
        
    end
end

return