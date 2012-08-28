function [numcomp alpha numcompT alphaT] = AFQ_MultiCompCorrection(data, method)
% Compute a multiple comparison correction for an AFQ data structure
%
%   [numcomp alpha numcompT alphaT] = AFQ_MultiCompCorrection(data)
%
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

switch(method)
    case 'simulate'
        % If the data is a single matrix of data then there is only one correction
        % to be done
        if size(data,1) > 2

                for ii = 1:100
                    x = randn(size(data,1),1);
                    [c(ii,:) p(ii,:)] = corr(x, data, 'rows', 'pairwise');
                end

                cMax = max(c');
                pMin = min(p');
                
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
                alphaT(ii) = 1 - (1 - 0.05)^(1 ./ numcompT(ii));
                
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
            alpha = 1 - (1 - 0.05)^(1 ./ numcomp);
        end
    case 'chevrud'
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
            alpha = 1 - (1 - 0.05)^(1 ./ numcomp);
            % Since data was not sent in separately for each tract we can assume
            % that the number of comparisons per tract is the same as the overall
            numcompT = numcomp; alphaT = alpha;
            
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
                alphaT(ii) = 1 - (1 - 0.05)^(1 ./ numcompT(ii));
                
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
            alpha = 1 - (1 - 0.05)^(1 ./ numcomp);
        end
end
return


