function params = CreateParamsStruct(c)
% Create a parameters structure from a cell array.

params = struct;
params.type = 'params version 1';
for ii = 1:2:length(c)
    params.(c{ii}) = c{ii+1};
end