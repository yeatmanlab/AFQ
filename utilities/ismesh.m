function L = ismesh(s)
% Check if a variable is an AFQ msh mesh structure
%
% L = ismsh(s)
%
% Returns a logical, true if it is a mesh and false if not

if exist('s','var') && ~isempty(s) && isstruct(s) && isfield(s,'tr') && isfield(s,'vertex')
    L = true;
else
    L = false;
end