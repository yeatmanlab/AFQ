function f = AFQ_TractProfileInterp(params,v,showmapping)

% function f = mapping(params,v,showmapping)
%
% <params> is 2 x N with the coordinates of the N anchor points
% <v> is a vector of data
%
% map <v> according to <params> and return the new vector.
%
% example:
% v = evalgaussian1d([20 5 1 0],1:100);
% v2 = mapping([.3 .7]',v);
% figure; hold on;
% plot(v,'r-');
% plot(v2,'b-');

if ~exist('showmapping','var') || isempty(showmapping)
    showmapping = false;
end

% calc
xval = linspace(0,1,length(v));


params = reshape(params,2,[]);  % first row is x-coordinates

% % enforce 0 and 1 bounds
% params = max(0,params);
% params = min(1,params);

% % sort according to x-values
% [params(1,:),ii] = sort(params(1,:));
% params(2,:) = params(2,ii);

% if not monotonically increasing, blow up
if ~all(diff(params(1,:)) > 0) || ~all(diff(params(2,:)) > 0)
    f = zeros(size(v));
    return;
end

% %
% if any(diff(params(1,:)) < .1) || any(diff(params(2,:)) < .1)
%   f = zeros(size(v));
%   return;
% end

% get the locations that we want
xvaldesired = interp1([0 params(2,:) 1],[0 params(1,:) 1],xval,'cubic');

% do the interpolation through the original vector to
% get the final vector
f = interp1(xval,v,xvaldesired,'cubic');

% Display mapping
if showmapping == 1
    figure(999); clf;
    plot(xval,xvaldesired,'ro-');
    xlabel('Anchor Point');
    ylabel('Mapping Destination');
    pause(.02);
end
