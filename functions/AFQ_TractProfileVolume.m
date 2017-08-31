function [Volume] = AFQ_TractProfileVolume(fg)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Calculate the volume along the tract profile.
%
% [Volume] = AFQ_TractProfileVolume(fg)
% or
% [TractProfile] = AFQ_TractProfileVolume(TractProfile)
% 
% If fg is a fiber group than volume will be calculated based on the number
% of unique coordinates at each node. If fg is a TractProfile structure
% than volume will be calculated based on the covariance matrix in the
% tract profile.

% Check if a tract profile or a fiber group was passed in
if exist('fg') && ~isempty(fg)
    if isfield(fg,'fibercov')
        method = 'cov';
    elseif isfield(fg,'fibers')
        method = 'fibers';
    end
else
    error('Pass in a fiber group or tract profile')
end

% Different calculation depending on the input
switch(method)
    case 'cov'
        % Check for empty fiber groups
        FGnotempty=zeros(length(fg),1);
        for jj = 1:length(fg)
            FGnotempty(jj) = ~isempty(fg(jj).fibercov);
        end
        FGnotempty = find(FGnotempty)';
        
        % Loop over the number of tract profiles
        for jj = FGnotempty
            % Extract the fiber covariance in matrix form
            fcov = fiberCovToMatrix(fg(jj));
            % Number of nodes for this Tract Profile
            numnodes = size(fcov,3);
            % Loop over nodes and calculate volume
            for ii = 1:numnodes
                S = fcov(:,:,ii);
                % Equation from -  http://en.wikipedia.org/wiki/Ellipsoid#Volume
                V(ii,jj) = (4/3)*pi*sqrt(det(S));
            end
            % Assign volume to the TractProfile
            fg(jj).fiberCovVolume = V(:,jj)';
        end
        % Assign the computed tract profile to the output variable
        Volume = fg;
    case 'fibers'
        
        % Check for empty fiber groups
        FGnotempty=zeros(length(fg),1);
        for jj = 1:length(fg)
            FGnotempty(jj) = ~isempty(fg(jj).fibers);
            numNodes(jj) = size(fg(jj).fibers{1},2);
            if numNodes(jj) ~= size(fg(jj).fibers{2},2);
                error('\nfiber group must be resampled to consistend number of nodes');
            end
        end
        FGnotempty = find(FGnotempty)';
        if length(fg)>1 && all(diff(numNodes))~=0
            error('\nfiber group must be resampled to consistend number of nodes');
        else
            numNodes=numNodes(1);
        end
        
        % Loop over the number of fiber groups
        for jj = FGnotempty
            % Concatenate fiber coordinates into a 3d array where each
            % entry in the 3rd dimension is a fiber
            f = cat(3,fg(jj).fibers{:});
            % Loop over each node and compute the number of unique
            % coordinates
            for ii = 1:numNodes
                Volume(ii,jj) = size(unique(round(squeeze(f(:,ii,:)))','rows'),1);
            end
        end
end