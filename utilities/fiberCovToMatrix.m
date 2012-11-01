function sigma = fiberCovToMatrix(fibercov)
% Convert the fiber covariance structure in a tract profile to matrix form
% The output is a 3x3xN matrix where N is the number of nodes in the fiber.
%
% Written by Jason D. Yeatman September 2012

if ~exist('fibercov','var') || isempty(fibercov)
    error('Please provide a fiber covariance structure');
elseif isstruct(fibercov) && isfield(fibercov,'fibercov')
    fibercov = fibercov.fibercov;
end
numnodes = size(fibercov,2);
% Put the entries in the right spots
for node = 1:numnodes
    sigma(:,:,node) = [fibercov(1:3, node)'; ...
        fibercov(2,node)', fibercov(4:5, node)';...
        fibercov(3,node)', fibercov(5,node)', fibercov(6, node)'];
end

return