function [alphaFWE, statFWE, clusterFWE, stats] = AFQ_MultiCompCorrection(data,y,alpha,method,cThresh)
% Compute a multiple comparison correction for Tract Profile data
%
%   [alphaFWE, statFWE, clusterFWE, stats] = AFQ_MultiCompCorrection(data,y,alpha, method,cThresh)
%
% There are 2 multiple comparison corrections implemented. Both account for
% the correlation structure in the data but in different ways.
%
% method = 'permutation' (default)
% This is an implementation of the permutation method described by Nichols
% and Holmes (2001). Nonparametric permutation tests for functional
% neuroimaging: A primer with examples. Human Brain Mapping.  This will
% return the faily wise error (FWE) corrected alpha value for pointwise
% comparisons.  It will also compute the FWE corrected cluster size at the
% user defined alpha.  This means that significant clusters of this size or
% greater are pass the multiple comparison threshold and do not need 
% further p-value adjustment.
%
% method = 'chevrud'
% This is an implementation of the multiple comparison correction proposed
% in Cheverud, J. M. (2001). A simple correction for multiple comparisons
% in interval mapping genome scans. Heredity (Edinb), 87(Pt 1), 52-58.  It
% calculates the number of independent variables in a dataset based on the
% varience in the eigenvalues of the correlation matrix for the data.  For
% Tract Diffusionp Profiles, values at nearby nodes are highly correlated
% and should not be treated as independent variables.  By taking into
% account the correlation between values at multiple points and on multiple
% tracts this alogorithm will determine a much more reasonable multiple
% comparison correction than the typical, overly conservative bonferroni
% correction
%
% Inputs:
% data     = Either a matrix of data for a single tract, or a matrix of data
%            for all the tracts combined.
% y        = A vector of either behavioral measurements or a binary 
%            grouping variable for which pointwise statistics will be
%            computed on the Tract Profile and the p-value adjusted for
%            mulltiple comparisons will be determined.  If y is a
%            continuous variable then correlations will be computed. If y 
%            is a binary vector then T-tests will be computed.
% alpha    = The desired alpha (pvalue) to adjust
% method   = 'permutation' or 'chevrud'. We strongly recomend 'permutation'
% cThresh  = For clusterwise corrections the threshold for computing a
%            cluster can be different than the desired alpha. For example
%            you can set a cluster threshold of 0.01 and then find clusters
%            that a large enough to pass FWE at a threshold of 0.05.
%
% Outputs:
% alphaFWE   = This is the alpha (p value) that corresponds after adjustment 
%              for multiple comparisons
% statFWE    = This is the value of the statistic corresponding to alphaFWE.
%              statFWE will either be a correlation coeficient or T-statistic
% clusterFWE = Clusters of points on a Tract Profile that are larger than
%              clusterFWE are significant at pvalue = alpha.
% stats      = A structure containing the results of each permutation
%
% Example:
%
% % Get a matrix of Tract FA Profile values for all the tracts
% data = AFQ_get(afq,'all vals', 'fa');
% % Compute the corrected p-value for ttests along the tract
% [alphaFWE statFWE] = AFQ_MultiCompCorrection(data,sub_group,0.05)
% 
% Written by Jason D. Yeatman, August 2012


if ~exist('method','var') || isempty(method)
    method = 'permutation';
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
% By default the cluster threshold is the same as alpha
if ~exist('cThresh') || isempty(cThresh)
    cThresh = alpha;
end
% If y is continues perform a correlation if binary perform a ttest
if ~exist('y','var') || isempty(y)
    y = randn(size(data,1),1);
    fprintf('No behavioral data provided so randn will be used');
    stattest = 'corr';
elseif length(y)==sum(y==0 | y==1) || length(y)==sum(y==1 | y==2)
    y = logical(y);
    stattest = 'ttest';
else
    stattest = 'corr';
end

switch(method)
    %% Permutation method
    case 'permutation'
        % number of permutations
        nperm = 1000;
        
        switch(stattest)
            case('corr')
                for ii = 1:nperm
                    % Shuffle the rows of the data
                    rows = Shuffle(1:size(data,1));
                    [stat(ii,:) p(ii,:)] = corr(y, data(rows,:), 'rows', 'pairwise');
                end
            case('ttest')
                for ii = 1:nperm
                    rows = Shuffle(y);
                    [~,p(ii,:),~,tstats] = ttest2(data(rows,:),data(~rows,:));
                    stat(ii,:) = tstats.tstat;
                end
        end
        % Sort the pvals and associated statistics such that the first
        % entry is the most significant
        stats.pMin = sort(min(p,[],2),'ascend');
        stats.statMax = sort(max(stat,[],2),'descend');
        % Find the corrected alpha and statistic
        alphaFWE = stats.pMin(round(alpha.*nperm));
        statFWE = stats.statMax(round(alpha.*nperm));
        
        %% Cluster threshold
        % If a cluster size is defined, also determine the significant
        % cluster size at the specified alpha value
        
        % Threshold the pvalue
        pThresh = p < cThresh;
        % Compute cluster size for each permutation
        for ii = 1:nperm
            % Find indices where significant clusters end.
            % The method used requires significant p-values to be included
            % between non-significant p-values. 0 are therefore added at 
            % both ends of the thresholded p-value vector 
            %(for cases when significant p-values are located at its ends)
            pThresh_ii=[0 pThresh(ii,:) 0];
            clusEnd = find(pThresh_ii == 0);
            % Compute the size of each cluster
            clusSiz = diff(clusEnd);
            % Find the maximum cluster size for permutation ii
            clusMax(ii) = max(clusSiz);
        end
        % Sort the clusters in descending order of significance
        stats.clusMax = sort(clusMax,'descend');
        % Find the corrected cluster size corresponding to alpha
        clusterFWE = stats.clusMax(round(alpha.*nperm));
        
    case 'chevrud'
        %% PCA method (Chevrud et al.)
        % If the data is a single matrix of data then there is only one correction
        % to be done
        if size(data,1) > 2
            
            % Compute the correlation matrix between the columns of the data
            if sum(isnan(data(:))) == 0
                c = corr(data);
            else
                c = corr(data,'rows','pairwise');
            end
            % Compute the eigenvalues of the correlation matrix
            s = eig(c);
            % Calculate the variance in the eigenvalues
            Vobs = var(s);
            % Count the number of colums in the data
            M = size(data,2);
            % This is the equation to calculate the number of independent
            % observations.
            %  See Chevrud (2000)
            numcomp = M*(1 - (M -1)*Vobs ./ (M^2));
            % This is the equation for the bonferroni correction with an alpha of
            % .05
            alphaFWE = 1 - (1 - 0.05)^(1 ./ numcomp);
            % Since data was not sent in separately for each tract we can assume
            % that the number of comparisons per tract is the same as the overall
            numcompT = numcomp; alphaTractFWE = alpha;
            
            % If a cell array of data for multiple tracts was sent in then we will
            % compute for each tract separately and then for the full dataset
        elseif isstruct(data)
            
            % Compute the number of independent comparisons for each tract
            % separately
            for ii = 1:length(data)
                % Compute the correlation matrix between the columns of the data
                c = corr(data(ii).FA,'rows','pairwise');
                % Compute the eigenvalues of the correlation matrix
                s = eig(c);
                % Calculate the variance in the eigenvalues
                Vobs = var(s);
                % Count the number of colums in the data
                M = size(data(ii).FA,2);
                % This is the equation to calculate the number of independent
                % observations.
                %  See Chevrud (2000)
                numcompT(ii) = M*(1 - (M -1)*Vobs ./ (M^2));
                % This is the equation for the bonferroni correction
                alphaTractFWE(ii) = 1 - (1 - 0.05)^(1 ./ numcompT(ii));
                
                clear c s Vobs M
            end
            % No do the correction for the full dataset rather than each fiber
            % group independently
            datall = horzcat(data(:).FA);
            % Compute the correlation matrix between the columns of the data
            c = corr(datall,'rows','pairwise');
            % Compute the eigenvalues of the correlation matrix
            s = eig(c);
            % Calculate the variance in the eigenvalues
            Vobs = var(s);
            % Count the number of colums in the data
            M = size(datall,2);
            % This is the equation to calculate the number of independent
            % observations.
            %  See Chevrud (2000)
            numcomp = M*(1 - (M -1)*Vobs ./ (M^2));
            % This is the equation for the bonferroni correction
            alphaFWE = 1 - (1 - 0.05)^(1 ./ numcomp);
        end
end

return


