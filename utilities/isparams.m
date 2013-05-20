function L = isparams(s)
% Check if a variable is a params structure

if isstruct(s) && isfield(s,'type') && strcmp('params version 1',s.type)
    L = true;
else
    L = false;
end