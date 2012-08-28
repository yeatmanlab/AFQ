function [alphaFWE statFWE clusterFWE stats] = AFQ_MultiCompCorrection(data,y,alpha, method, cluster)
% Compute a multiple comparison correction for Tract Profile data
%
%   [alphaFWE statFWE clusterFWE] = AFQ_MultiCompCorrection(data,y,alpha, method)
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
% greater are significant at alpha with p-value adjustment.
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
% data     = Either a matrix of data for a single tract, a matrix of data
%            for all the tracts combined, or a cell array of data where
%            each cell is the data for a single tract.
% y        =
% alpha    =
% method   =
%
% Outputs:
% numcomp = This is the total number of independent comparisons for the
%           full dataset.  If the data is a cell array of data for
%           different fiber groups then numcomp represents the number of
%           comparisons for all the data concatinated together
% alpha   = This is the alpha (p value) that corresponds to a alpha of 0.05
%           after the correction for the number of independent comparisons.
%           If the data is a cell array of data for different fiber groups
%           then alpha represents the proper alpha if comparisons are done
%           for every location on every tract
% numcompT= The same as numcomp but done independently for each tract.
%           This is apropriate to use if statistics are only being done
%           for a single tract.
% alphaT  = This is the same as alpha but done independently for each tract

if ~exist('method','var') || isempty(method)
    method = 'permutation';
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
if ~exist('cluster','var') || isempty(cluster)
    cluster = 1;
end
% If y is continues perform a correlation if binary perform a ttest
if ~exist('y','var') || isempty(y)
    y = randn(size(data,1),1);
    fprintf('No behavioral data provided so randn will be used');
    stat = 'corr';
elseif length(y)==sum(y==0 | y==1) || length(y)==sum(y==1 | y==2)
    y = logical(y)
    stat = 'ttest'
end

switch(method)
    %% Permutation method
    case 'permuatation'
        % number of permutations
        nperm = 1000;
        
        switch(stat)
            case('corr')
                for ii = 1:nperm
                    % Shuffle the rows of the data
                    rows = Shuffle(1:size(data,1));
                    [stat(ii,:) p(ii,:)] = corr(y, data(rows,:), 'rows', 'pairwise');
                end
            case('ttest')
                for ii = 1:nperm
                    rows = Shuffle(y);
                    for jj = 1:size(data,2)
                        [~,p(ii,jj),~,tstats] = ttest2(data(rows,jj),data(~rows,jj));
                        stat(ii,jj) = tstats.tstat;
                    end
                end
        end
        % Sort the pvals and associated statistics such that the first
        % entry is the most significant
        stat.pMin = sort(min(p,[],2),'ascend');
        stat.statMax = sort(max(stat,[],2),'descend');
        % Find the corrected alpha and statistic
        alphaFWE = stat.pMin(round(alpha.*nperm));
        statFWE = stat.statMax(round(alpha.*nperm));
        
        %% Cluster threshold
        % If a cluster size is defined, also determine the significant
        % cluster size at the specified alpha value
        
        % Threshold the pvalue
        pThresh = p < alpha;
        % Compute cluster size for each permutation
        for ii = 1:nperm
            % Find indices where significant clusters end
            clusEnd = find(pThresh(ii,:) == 0);
            % Compute the size of each cluster
            clusSiz = diff(clusEnd);
            % Find the maximum cluster size for permutation ii
            clusMax(ii) = max(clusSiz);
        end
        % Sort the clusters in descending order of significance
        stat.clusMax = sort(clusMax,'descend');
        % Find the corrected cluster size corresponding to alpha
        clusterFWE = clusMax(round(alpha.*nperm))
        
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


