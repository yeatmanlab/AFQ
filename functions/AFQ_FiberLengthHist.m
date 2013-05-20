function [Lnorm Lmm]=AFQ_FiberLengthHist(fg, showHist)
% Calculate the distribution of fiber lengths for a fiber group
%
% [Lnorm Lmm]=AFQ_FiberLengthHist(fg, [showHist = fales])
%
% Inputs:
% fg       = The input fiber group
% showHist = 0 Don't display a histagram of lengths, 1 do display histogram
%
% Outputs:
% Lnorm    = Normalized length of each fiber, in units of z-score.
% Lmm      = The length of each fiber in mm.
%
% Written by Jason D. Yeatman 11/5/2011

if ~exist('fg','var') || isempty(fg) || ~isfield(fg,'fibers') || isempty(fg.fibers)
    fprintf('\nNo fibers in fiber group')
    Lnorm=[]; Lmm=[];
    return
end
if ~exist('showHist','var') || isempty(showHist)
    showHist = 0; %Default to not showing a figure
end
% Calculate the length in mm of each fiber in the fiber group
Lmm=cellfun('length',fg.fibers);
% Calculate distribution of lengths
[Lnorm, Mu, Sigma]= zscore(Lmm);
% Show a histagram
if showHist==1
    figure;hold on;
    [n xout]=hist(Lmm,20);
    bar(xout,n,'FaceColor','b','EdgeColor','k');
    axis([min(xout) max(xout) 0 max(n)]);
    xlabel('Fiber Length (mm)');
    plot([Mu Mu],[0 max(n)],'r','linewidth',2);
    plot([Mu-Sigma Mu-Sigma],[0 max(n)],'--r');
    plot([Mu+Sigma Mu+Sigma],[0 max(n)],'--r');
    plot([Mu-2*Sigma Mu-2*Sigma],[0 max(n)],'--r');
    plot([Mu+2*Sigma Mu+2*Sigma],[0 max(n)],'--r');
end

return