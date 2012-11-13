function [TractProfile V] = AFQ_TractProfileVolume(TractProfile)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Calculate the volume of the ellipse defined by the covariance matrix at
% each node of the tract profile
%
% [TractProfile V] = AFQ_TractProfileVolume(TractProfile)

% Check for empty fiber groups
FGnotempty=zeros(length(TractProfile),1);
for jj = 1:length(TractProfile)
    FGnotempty(jj) = ~isempty(TractProfile(jj).fibercov);
end
FGnotempty = find(FGnotempty)';

% Loop over the number of tract profiles
for jj = FGnotempty
    % Extract the fiber covariance in matrix form
    fcov = fiberCovToMatrix(TractProfile(jj));
    % Number of nodes for this Tract Profile
    numnodes = size(fcov,3);
    % Loop over nodes and calculate volume
    for ii = 1:numnodes
        S = fcov(:,:,ii);
        % Equation from -  http://en.wikipedia.org/wiki/Ellipsoid#Volume
        V(ii,jj) = (4/3)*pi*sqrt(det(S));
    end
    % Assign volume to the TractProfile
    TractProfile(jj).fiberVolume = V(:,jj)';
end