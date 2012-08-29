function [coeff, score, latent] = AFQ_pca(afq, valname, group)
% Perform principal components analysis on AFQ Tract Profiles
%
% [coeff, score, latent] = AFQ_pca(afq, valname, group)
%
% See help princomp for more information on coeff, score and latent.
%
% Input:
% afq     - afq structure.  See AFQ_Create
% valname - Name of the value to perform pca on. Default valname = 'fa'
% group   - Compute PCA on values from either 'patients' or 'controls'. If
%           group is empty then PCA will be computed on all the values
%
% Example:
% 
% % Perform PCA on the controls' Tract FA Profiles
% [coeff, score, latent] = AFQ_pca(afq, 'fa', 'controls');
%
% Written by Jason D. Yeatman, August, 2012


if ~exist('afq','var') || isempty(afq) || ~isafq(afq)  
    error('Please provide and afq structure')
end
if ~exist('valname','var') || isempty(valname)
    valname = 'fa';
    fprintf('Performing PCA on FA values')
end

% Concatinate all the data into an MxN matrix where each row is a
% subject and each column is a point on a tract profile. If there are 10
% subjects with 20 tracts and 100 points per tract than data will be a
% 10x2000 matrix
if ~exist('group','var') || isempty(group)
    data = AFQ_get(afq, 'vals', valname);
else
    data = AFQ_get(afq, 'vals', valname, group);
end

% Perform PCA
if sum(isnan(data(:))) > 0
    % identify columns with nans
    nancol = sum(isnan(data) > 0);
    % remove columns with nans
    fprintf('\nremoving %.0f columns  because they contain nans\n',sum(nancol));
    dataNoNan = data(:,~nancol);
    % PCA on all columns with no nans
    [coeff, score, latent] = princomp(dataNoNan);
else
    [coeff, score, latent] = princomp(data);
end