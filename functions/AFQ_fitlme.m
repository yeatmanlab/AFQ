function lme = AFQ_fitlme(afq, response, predictors, fgname, m, c, outliers)
% Fit mixed linear model to data in afq structure
%
% lme = AFQ_fitlme(afq, response, predictors, fgname, c, outliers)
%
% Inputs:
%
% afq        - afq structure. See AFQ_run.m
% response   - property to use as response variable (y). default 'fa'
% predictors - names of variables that should be predictors. These names
%              should correspond to variables in afq.metadata
% fgname     - name of the fiber group
% m          - binary. predict the tract mean or each node along the tract
% c          - binary. should predictors be defined as categorical
% outliers   - binary. should outliers be excluded? Will look in
%              afq.metadata.outiers
%
% fgnames = AFQ_get(afq, 'fgnames')
% [outliers ,afq] = AFQ_outliers(afq, {'dki_MK', 'dki_MD'}, 4, 40)
% for ii = 1:length(fgnames)
%    results{ii} = AFQ_fitlme(afq, 'dki_MD', 'session', fgnames{ii}, 1, 1)
% end

if ~exist('response', 'var') || isempty(response)
    response = 'fa';
end
if ~exist('predictors', 'var') || isempty(predictors)
    predictors = {'session'};
    fprintf('\nSince no predictors were defined we will use session\n')
elseif ischar(predictors)
    tmp = predictors; clear predictors;
    predictors{1} = tmp;
end

if ~exist('outliers', 'var') || isempty(outliers) || outliers == 0
    exclude = zeros(length(afq.sub_dirs),1);
else
    exclude = AFQ_get(afq, 'metadata', 'outliers');
end

% Get the response values
y = AFQ_get(afq, fgname, response);

% Take the mean if interested in analyzing the tract average (else, loop
% through tract nodes - TO DO, add option to return mc corrected pvals)
if ~exist('m', 'var') || m == 1
    y = nanmean(y,2);
end

% Check for nans
na = isnan(y);
if sum(na,2) > 0
    fprintf('Removing subject %s due to nans\n',afq.sub_names{na});
    exclude = exclude | na;
end

% Get the predictors
for ii = 1:length(predictors)
    x(:,ii) = AFQ_get(afq, 'metadata', predictors{ii});
end

% Remove outliers if desired
y = y(~exclude, :);
x = x(~exclude, :);

% Loop over nodes, or analyze a single mean value if m == 1
for n = 1:size(y,2)
    yy = y(:,n);
    
    % Create a table with the data
    cnames = predictors; cnames = horzcat({'subject'},{response},predictors)
    if c == 1
        d = table(afq.sub_names(~exclude,:), yy, categorical(x),...
            'VariableNames', cnames);
    else
        d = table(afq.sub_names(~exclude,:), yy, x,...
            'VariableNames', cnames);
    end
    
    % Concatenate a string with all predictor names
    p = [];
    for ii = 1:length(predictors)
        p = strcat(p, sprintf('%s + ',predictors{ii}));
    end
    % Fit the model
    model = sprintf('%s ~ %s (%s | subject)', response, p, p(1:end-1))
    %model = sprintf('%s ~ %s + (1 | subject)', response, predictors{:})
    
    lme{n} = fitlme(d, model);
    
end