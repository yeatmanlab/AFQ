function [curvWA torsWA TractProfile] = AFQ_ParamaterizeTractShape(fg, TractProfile)
% Paramaterize curvature and torsion along a tract profile
%
% [curvWA torsWA TractProfile] = AFQ_ParamaterizeTractShape(fg, TractProfile)
%
% Inputs:
% fg           - Fiber Group. This can be a structure containing multiple
%                fiber groups where each fiber group is an entry in the
%                structure
% TractProfile - Tract Profile for the fiber group. If one is not entered
%                then it will be calculated
%
% Outputs:
% curvWA       - Curvature of the tract core calculated by taking a
%                weighted average of curvature on each fiber.
% torsWA       - Torsion of the tract core calculated by taking a weighted
%                average of torsion on each fiber.
% TractProfile - Tract Profile containing these measurements
%
% Written by Jason D. Yeatman September 2012

%% Argument checking
if ~exist('fg','var') || isempty(fg)
    error('Please provide a fiber group')
elseif isfield(fg,'subgroup') && isfield(fg,'subgroupNames') && length(fg.subgroupNames) > 1
    fg = fg2Array(fg);
end

if ~exist('TractProfile','var') || isempty(TractProfile)
    % Build a tract profile of the fiber group if one was not passed in
    numNodes = 100;
    % If there are multiple fiber groups stored within fg then loop over them
    for ii = 1:length(fg)
        [SuperFiber, fg(ii)]= dtiComputeSuperFiberRepresentation(fg(ii), [], numNodes);
        TractProfile(ii) = AFQ_CreateTractProfile('name',fg.name,'superfiber',SuperFiber);
    end
elseif length(TractProfile) == length(fg)
    for ii = 1:length(fg)
        % Number of nodes in the tract profile
        numNodes(ii) = length(TractProfile(ii).coords.acpc);
        % Resample and reorient fibers so that each fiber starts and ends in the
        % same place
        if ~isempty(fg(ii).fibers) && numNodes(ii) > 0
            fg(ii) = dtiReorientFibers(fg(ii), numNodes(ii));
        end
    end
    % Sometimes a group has 0 nodes because it is empty. This will mess up
    % future computations so we se numNodes to the maximum number of nodes.
    % Note that all groups probably should have the same number of nodes
    numNodes = max(numNodes);
elseif length(TractProfile) ~= length(fg)
    error('TractProfile and fg must have the same number of fiber groups')
end

%% Loop over the fiber groups in fg and calculate torsion and curvature
% Allocate data matrices
curvWA = nan(numNodes,length(fg));
torsWA = nan(numNodes,length(fg));
% Check for empty fiber groups
FGnotempty=zeros(length(fg),1);
for jj = 1:length(fg)
    FGnotempty(jj) = ~isempty(fg(jj).fibers) && ~isempty(TractProfile(jj).fibercov);
end
FGnotempty = find(FGnotempty)';
for jj = FGnotempty
    % Initialize variables
    numfibers = length(fg(jj).fibers);
    curv = zeros(numNodes,numfibers);
    tors = zeros(numNodes,numfibers);
    weights = zeros(numNodes,numfibers);
    
    % Calculate curvature and torsion for each node on each fiber
    for ii = 1:numfibers
        [~,~,~,curv(:,ii),tors(:,ii)] = frenet2(fg(jj).fibers{ii}(1,:)',fg(jj).fibers{ii}(2,:)',fg(jj).fibers{ii}(3,:)');
    end
    
    % Each fiber is represented by numberOfNodes, so can easily loop over 1st
    % node, 2nd, etc...
    fc = horzcat(fg(jj).fibers{:})';
    
    % Calculate the weights to be applied to each fibers curvature estimates.
    % These weights are based on the fiber's mahalanobis distance from the
    % fiber core. They will be used to weight that fibers contribution to the
    % core measurement
    for node=1:numNodes
        % Get the covariance structure from the tract profile
        fcov = TractProfile(jj).fibercov;
        % Compute gaussian weights y = mvnpdf(X,mu,SIGMA);
        % Returns the density of the multivariate normal distribution with zero
        % mean and identity covariance matrix, evaluated at each row of X.
        X=fc((1:numNodes:numfibers*numNodes)+(node-1), :);
        sigma = [fcov(1:3, node)'; ...
            0 fcov(4:5, node)';...
            0 0 fcov(6, node)'];
        sigma = sigma + sigma' - diag(diag(sigma));
        mu    = TractProfile(jj).coords.acpc(:, node)';
        % Weights for the given node.
        weights(node, :) = mvnpdf(X,mu,sigma)';
    end
    
    % The weights for each node should add up to one across fibers.
    weightsNormalized = weights./(repmat(sum(weights, 2), [1 numfibers]));
    
    % Weight each curvature and torsion calculation by its distance from the
    % core
    curvWA(:,jj) = sum(curv.*weightsNormalized,2);
    torsWA(:,jj) = sum(tors.*weightsNormalized,2);
    
    % Assign values to the tract profile structure
    TractProfile(jj).fiberCurvature = curvWA(:,jj);
    TractProfile(jj).fiberTorsion   = torsWA(:,jj);
    
end