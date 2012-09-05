function [params,mapv1,mapv1b] = AFQ_TractProfileAlign(v1,v2,numanchor)
% Align two Tract Profiles by stretching one to match the other
%
% function [params,mapv1,mapv1b] = AFQ_TractProfileAlign(v1,v2,numanchor)
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
% numanchor - The number of anchor points on v1 to map to v2. 
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

% define options for the optimization
options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun',1e-6,'TolX',1e-6);

% define seeds for the parameter search
params0 = linspace(0,1,2+numanchor);
params0 = repmat(params0(2:end-1),[2 1]);

% figure out seed for scale and offset with regression
h = pinv([v1' ones(length(v1),1)])*v2';  % first parameter is scale, second is offset

% construct final seed
params0 = [h(:)' params0(:)'];

% define bounds paramslb = zeros(1,numanchor); paramsub =
% ones(1,numanchor);
paramslb = [-Inf -Inf zeros(1,2*numanchor)];
paramsub = [Inf Inf ones(1,2*numanchor)];

% Fit parameters OLD: sqrt(1-nanreplace(corr(mapping(a,v1)',v2'),-1))
[params,d,d,exitflag,output] = ...
    lsqnonlin(@(params) v2 - (params(1)*AFQ_TractProfileInterp(params(3:end),v1) + params(2)) , ...
    params0,paramslb,paramsub,options);
assert(exitflag > 0);  % ensure good convergence

% Apply stretching to v1 so it matches v2
mapv1 = AFQ_TractProfileInterp(params(3:end),v1);
% Apply scale and offset
mapv1b = params(1)*mapv1 + params(2);

