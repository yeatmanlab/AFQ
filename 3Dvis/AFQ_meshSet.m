function msh = AFQ_meshSet(msh, param, varargin)
% Set fields in the AFQ msh structure
%
%

switch(param)
    case('vertices')
        if nargin < 3, fprintf('Please supply vertex name'); end
        if isfield(msh.vertex, varargin{1})
           msh.vertices = msh.vertex.(varargin{1});
        end  
end