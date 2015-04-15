function rgbvals = lineargradient(rgbends,nsamples)
%
% example:
% rgbends = [.7 .6 .6; .7 0 0]; nsamples = 256;
% rgbvals = lineargradient(rgbends,nsamples)

rgbvals = vertcat(linspace(rgbends(1,1),rgbends(2,1),nsamples),...
    linspace(rgbends(1,2),rgbends(2,2),nsamples),...
    linspace(rgbends(1,3),rgbends(2,3),nsamples))';