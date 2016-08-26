function lme = AFQ_fitlme(afq, response, predictors, fgname, c, outliers)
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
% c          - binary. should predictors be defined as categorical
% outliers   - binary. should outliers be excluded? Will look in
%              afq.metadata.outiers
%
% fgnames = AFQ_get(afq, 'fgnames')
% for ii = 1:length(fgnames)
%    lme(ii) = AFQ_fitlme(afq, 'dki_FA', 'session', fgnames{ii}, 1, 1)
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

% Get the predictors
for ii = 1:length(predictors)
    x(:,ii) = AFQ_get(afq, 'metadata', predictors{ii});
end

% Remove outliers if desired
y = y(~exclude, :);
x = x(~exclude, :);

% Create a table with the data
if c ==1
    d = table(afq.sub_names(~exclude,:), nanmean(y,2), categorical(x), 'VariableNames', {'subject', response, predictors{:}});
else
    d = table(afq.sub_names(~exclude,:), nanmean(y,2), categorical(x), 'VariableNames', {'subject', response, predictors{:}});
end
% Fit the model
model = sprintf('%s ~ %s + (%s | subject)', response, predictors{:}, predictors{:})
%model = sprintf('%s ~ %s + (1 | subject)', response, predictors{:})

lme = fitlme(d, model)