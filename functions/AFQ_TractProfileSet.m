function TractProfile = AFQ_TractProfileSet(TractProfile,param,varargin)
% Set fields of the afq TractProfile structure

switch(param)
    % Set coordinates
    case{'coords' 'coordsacpc' 'acpccoords'}
        TractProfile.coords.acpc = varargin{1};
    case{'vals'}
        % Name of the value is the first argument
        valname = varargin{1};
        % The values are the second
        vals = varargin{2};
        % Set the values
        TractProfile.vals.(valname) = vals;
        
end