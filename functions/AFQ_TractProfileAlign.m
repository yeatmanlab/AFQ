function [params,mapv1,mapv1b,err] = AFQ_TractProfileAlign(v1,v2,numanchor,errmetric,show)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Align two Tract Profiles by stretching one to match the other
%
% function [params,mapv1,mapv1b,err] = AFQ_TractProfileAlign(v1,v2,numanchor)
%
% This function will coregister the Tract Profile in vector 1 (v1) to the
% Tract Profile in vector 2 (v2).  This is done by using searching for an
% anchor point on v1 and a mapping of this anchor point to an equivalent
% point on v2. v1 is then interpolated with respect to this anchor point so
% that the two vectors match.  The search may have problems with local
% minimums when the number of anchor points is > 1.
%
% Inputs:
%
% v1        - An input vector to coregister
% v2        - Target vector. v1 will be coregisterd to v2.
% numanchor - The number of anchor points on v1 to map to v2. There may be
%             a problem with local minimum if you use more than 1
% errmetric - Error metric for fitting. Either 'R2' or 'MSE'
% show      - Plot out fitting proceedure. [True/False]
%
% Outputs:
%
% params    - Parameters of the transform.  These can be plugged into
%             AFQ_TractProfileInterp to resample v1.
% mapv1     - v1 coregistered to v2 just based on stretching the vector
%             along the x axis
% mapv1b    - v1 coregistered to v2 based on stretching it scaling it and
%             offsetting it
%
% Example:
%
% v1 = evalgaussian1d([20 5 1 0],1:100); v2 = evalgaussian1d([40
% 20 3 0],1:100); [params,mapv1,mapv1b] = AFQ_TractProfileAlign(v1,v2,2);
% figure; hold on; plot(v1,'r-'); plot(v2,'g-'); plot(mapv1,'b-');
% plot(mapv1b,'c-');

if ismatrix(v1) && (~exist('v2','var') || isempty(v2))
    v2 = nanmean(v1);
    nsub = size(v1,1);
else
    nsub = 1;
end

if ~exist('numanchor','var') || isempty(numanchor)
    numanchor = 1;
end

if ~exist('show','var') || isempty(show)
    show = 0;
end
if ~exist('errmetric','var') || isempty(errmetric)
    errmetric = 'mse';
end
% define options for the optimization
options = optimset('FunValCheck','on','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun',1e-6,'TolX',1e-6);

% define bounds paramslb = zeros(1,numanchor); paramsub =
% ones(1,numanchor);
paramslb = [-Inf -Inf zeros(1,2*numanchor)];
paramsub = [Inf Inf ones(1,2*numanchor)];

for ii = 1:nsub
    if sum(isnan(v1(ii,:))) > 0
        err(ii) = 1;
        mapv1(ii,:) = v1(ii,:);
        mapv1b(ii,:) = v1(ii,:);
        continue
    else
        err(ii) = 0;
        
    end
    % define seeds for the parameter search
    params0 = linspace(0,1,2+numanchor);
    params0 = repmat(params0(2:end-1),[2 1]);
    
    % figure out seed for scale and offset with regression
    % The first parameter is scale, second is offset
    h = pinv([v1(ii,:)' ones(length(v1(ii,:)),1)])*v2';
    
    % construct final seed
    params0 = [h(:)' params0(:)'];
    
    % Fit parameters OLD: sqrt(1-nanreplace(corr(mapping(a,v1(ii,:))',v2'),-1))
    [params,d,d,exitflag,output] = ...
        lsqnonlin(@(params) costfun(v1(ii,:),v2,params,errmetric,show) , ...
        params0,paramslb,paramsub,options);
    % ensure good convergence
    assert(exitflag > 0);
    
    % Apply stretching to v1 so it matches v2
    mapv1(ii,:) = AFQ_TractProfileInterp(params(3:end),v1(ii,:));
    % Apply scale and offset
    mapv1b(ii,:) = params(1)*mapv1(ii,:) + params(2);
end

% Plot the origional vector and the transformed
if show == 1
    figure;hold;
    plot(v2,'k');plot(v1,'b');plot(mapv1,'r');
    legend({'Destination' 'Origional' 'Transformed'})
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = costfun(v1,v2,params,metric,show,L)
% Cost function for fitting algorithm. Params(1) is the gain factor.
% Params(2) is the offset. Params(3:end) are anchor points and mappings.

if ~exist('metric','var') || isempty(metric)
    metric = 'mse';
end
if ~exist('L','var') || isempty(L)
    L = 0;
end
% Use coeficient of determination (R2) as error metric
switch(metric)
    case{'R2' 'r2'}
        v1Out = (params(1)*AFQ_TractProfileInterp(params(3:end),v1, show) + params(2));
        R2 = sum((v2-v1Out).^2)./sum((v2 - mean(v2)).^2);
        % Error is defined as the percent of unexplained variance. We take the sqrt
        % because lsq nonline squares error by default
        err = sqrt(1 - R2);
        
    case{'mse','MSE'}
        % Squared difference between input and output
        err = v2 - (params(1)*AFQ_TractProfileInterp(params(3:end),v1, show) + params(2));
        % Calculate the difference between the anchor point and the
        % mapping
        d = params(3:2:end) - params(4:2:end);
        % Create a second error metric that is the sum of the differences
        % scaled by some weight (L).
        err2 = sum(d.^2)*L;
        % Add this to the mean squared error
        err = abs(err)+err2;
       
end

