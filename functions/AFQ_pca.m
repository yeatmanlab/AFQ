function [coeff, score, subMeans, latent, tsquared, R2] = AFQ_pca(afq, valname, group, demean, fgnums)
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
% demean  - Whether or not to remove subject means before PCA
% fgnums  - Which tracts to analyze
% Output:
% coeff    - principal component coefficients, also known as loadings.Rows
%            of X correspond to observations, columns to variables.
%            coeff is a P-by-P matrix, each column containing coefficients
%            for one principal component.  The columns are in order of
%            decreasing component variance.
% score    - principal component scores, i.e., the representation of X in 
%            the principal component space.  Rows of SCORE correspond to 
%            observations, columns to components.
% latent   - returns the principal component variances, i.e., the 
%            eigenvalues of the covariance matrix. The variance explained
%            by each component can be calculated with:
%            R2 = cumsum(latent)./sum(latent);
% subMeans - Subject means.
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
if ~exist('demean','var') || isempty(demean)
    demean = 0;
end
if ~exist('fgnums', 'var') || isempty(fgnums)
    fgnums = 1:length(AFQ_get(afq, 'fgnames'));
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
% remove unwanted tracts
fgnodes = [];
for ii = 1:length(fgnums)
    fgnodes = horzcat(fgnodes,fgnums(ii)*100-99:fgnums(ii)*100);
end
data = data(:,fgnodes);
% Compute each subjects mean value
subMeans = nanmean(data,2);
% Remove each subjects mean
if demean == 1
    data = bsxfun(@minus,data,subMeans);
end

% Compute tract means
fgnames = AFQ_get(afq, 'fgnames');
fgnames = fgnames(fgnums)
nodes = 21:80;
for ii = 1:length(fgnames)
    vals = AFQ_get(afq, fgnames{ii}, valname);
    tmeans(:,ii) = nanmean(vals(:,nodes),2);
end

% Perform PCA
if sum(isnan(data(:))) > 0
    % identify columns with nans
    nancol = sum(isnan(data) > 0);
    % remove columns with nans
    fprintf('\nremoving %.0f columns  because they contain nans\n',sum(nancol));
    dataNoNan = data(:,~nancol);
    % PCA on all columns with no nans
    [coeff, score, latent, tsquared, R2] = pca(dataNoNan);
else
    [coeff, score, latent, tsquared, R2] = pca(data);
    [tcoeff, tscore, tlatent, ttsquared, tR2] = pca(tmeans);
end