function [SuperFiber, weights] = AFQ_FiberTractGaussian(fg, numberOfNodes, M)
% Represent fiber group as a 3d gaussian
%
% [SuperFiber, weights] = AFQ_FiberTractGaussian(fg,[numberOfNodes=30],[M='mean'])
% 
% Inputs:
%        fg           - fiber group structure
%        numberofNodes- sample the fiber group into this many points
%        M            - represent central tendency with 'mean' or 'median'
% Outputs:
%       Superfiber  - a structure describing the mean (core) trajectory and
%                     the spatial dispersion in the coordinates of fibers
%                     within a fiber group (3d covariance matrix).
%
%       weights  -    numberOfNodes by numberOfFibers array of weights
%                     denoting, the mahalanobis distance of each node in 
%                     each fiber from the fiber tract core
%  
%  Example:
%
%  AFQdata = '/home/jyeatman/matlab/svn/vistadata/AFQ'
%  fg = dtiReadFibers(fullfile(AFQdata,'fibers','Arcuate_L.mat'));
%  numberOfNodes = 30; M = 'mean';
%  [SuperFiber, weights] = AFQ_FiberTractGaussian(fg, numberOfNodes, M)
% 
%  Written by Jason D. Yeatman pm 12/8/2011.  
%  Based on code fromdtiFiberGroupPropertyWeightedAverage
% 

if notDefined('fg'), error('Fiber group required'); end
if notDefined('numberOfNodes'), numberOfNodes=100; end
if notDefined('M'), M = 'mean'; end
% Number of fibers in the whole fiber group
numfibers = size(fg.fibers, 1);
% This function will resample the fibers to numberOfNodes and will also
% reorient some fibers, so the notion of "first" and "last" may end up
% converted. [fg] returned in the line below will be resampled and reoriented.  
[SuperFiber, fg] = dtiComputeSuperFiberRepresentation(fg, [], numberOfNodes, M);
% Each fiber is represented by numberOfNodes, so can easily loop over 1st
% node, 2nd, etc...
fc = horzcat(fg.fibers{:})'; 
% Preallocate weights when you understand its size
% weights = zeros(numberOfNodes,???)
% Compute weights
weights=zeros(numberOfNodes,numfibers);
for node=1:numberOfNodes
    % Compute gaussian weights y = mvnpdf(X,mu,SIGMA);
    % Returns the density of the multivariate normal distribution with zero
    % mean and identity covariance matrix, evaluated at each row of X.
    X=fc((1:numberOfNodes:numfibers*numberOfNodes)+(node-1), :);
    sigma = [SuperFiber.fibervarcovs{1}(1:3, node)'; ...
        0 SuperFiber.fibervarcovs{1}(4:5, node)';...
        0 0 SuperFiber.fibervarcovs{1}(6, node)'];
    if rank(sigma)<3
        continue
    end
    sigma = sigma + sigma' - diag(diag(sigma));
    mu    = SuperFiber.fibers{1}(:, node)';
    % Weights for the given node.Are calculated based on Mahalanobis
    % distance
    d=bsxfun(@minus,X,mu); % each point minus the mean
    % This calculate tha mahalanobis distance of each point on each fiber
    % from the tract core
    weights(node,:) = sqrt(dot(d/(sigma), d,2))';
end

return
