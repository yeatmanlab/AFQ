function [fa, md, rd, ad, cl, SuperFiber, fgClipped, cp, cs, fgResampled] = ...
    cp_dtiComputeDiffusionPropertiesAlongFG(fg, dt, roi1, roi2, numberOfNodes, dFallOff)

%% This function computes a weighted average of a diffusion metric (FA/MD/RD/AD) in a segment of the cerebellar peduncles. 
%  It is based on the original "dtiComputeDiffusionPropertiesAlongFG" (vistasoft) -- no modifications were made.

%  From a fiber group (fg), and diffusion data (dt), compute the weighted value of a diffusion metric (taken from dt)
%  along a segment of the fiber group -- between two ROIs sampled at a specified number of nodes.

%  Original vistasoft function written by ER - December 2009
%  Adapted by Maya Yablonski - May 2016 
%  Edited for publication - January 2019

%  Copyright Stanford Vista Team 2008

%% Input variables
%  fg            - Fiber group structure, here: the cerebellar peduncles.
%  dt            - dt6.mat file or a nifti image.  
%                  If a nifti image is passed in then only 1 value will be 
%                  output and others will be empty
%  roi1          - First ROI for the fg
%  roi2          - Second ROI for the fg
%  numberOfNodes - Number of samples taken along each fg
%  dFallOff      - Rate of fall off in weight with distance. 


%% Argument checking

if notDefined('dFallOff'), dFallOff = 1; end

%  Check if the input is a nifti image rather than a dt6 file.
if isfield(dt,'qto_ijk')
    valname = 'image';
else
    valname = 'famdadrdShape';
end

%% Compute weighted average of a diffusion metric along the fiber group
%  If two ROIs are passed in clip the fiber group to the portion that spans between the ROIs
if ~notDefined('roi1') && ~notDefined('roi2')
    fgClipped = cp_dtiClipFiberGroupToROIs(fg,roi1,roi2);
    % Compute weighted averages for eigenvalues along clipped fiber tract
    [myValsFgWa, SuperFiber, ~, ~, fgResampled] = ...
        dtiFiberGroupPropertyWeightedAverage(fgClipped, dt, numberOfNodes, valname, dFallOff);
else
    % Compute weighted averages for eigenvalues along full fiber tract
    [myValsFgWa, SuperFiber, ~, ~, fgResampled] = ...
        dtiFiberGroupPropertyWeightedAverage(fg, dt, numberOfNodes, valname, dFallOff);
    % There is no clipped fiber group
    fgClipped = nan;
end

%% Pull out specific diffusion metrics
if strcmp(valname,'famdadrdShape')
    fa = myValsFgWa(:, 1);
    md = myValsFgWa(:, 2);
    ad = myValsFgWa(:, 3);
    rd = myValsFgWa(:, 4);
    cl = myValsFgWa(:, 5);
    cp = myValsFgWa(:, 6);
    cs = myValsFgWa(:, 7);
elseif strcmp(valname,'image')
    % If a nifti image was put in then just put the image values into the FA
    % variable and leave the other variables empty
    fa = myValsFgWa(:, 1);
    md = nan(numberOfNodes,1);
    ad = nan(numberOfNodes,1);
    rd = nan(numberOfNodes,1);
    cl = nan(numberOfNodes,1);
    cp = nan(numberOfNodes,1);
    cs = nan(numberOfNodes,1);
end

return