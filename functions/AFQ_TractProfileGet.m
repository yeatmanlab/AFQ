function val = AFQ_TractProfileGet(TractProfile,param,varargin)
% Get parameters from an afq tract profile
%
% AFQ_TractProfileGet(TractProfile,param,varargin)
%
%

% remove spaces and upper case
param = mrvParamFormat(param);

switch(param)
    % Get 
    case{'vals'}
        % The name of the value to be returned
        if nargin == 2
            valname = 'fa';
        else
            valname = varargin{1};
        end
        % Get the values
        if isfield(TractProfile.vals,valname)
            val = TractProfile.vals.(valname);
        else
            error('Tract Profile value %s does not exist',valname);
        end
    case{'coords' 'coordsacpc' 'acpc'}
        val = TractProfile.coords.acpc;
        
end