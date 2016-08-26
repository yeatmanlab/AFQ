function [outliers, afq] = AFQ_outliers(afq, properties, thresh, num)
% Detect outliers who have tracts with unexpected values
%
% outliers = AFQ_outliers(afq, properties, thresh, num)
%
% Inputs
%
% afq        - afq structure. See AFQ_run
% properties - properties upon which to detect outliers. Cell array
% thresh     - distance (z-score) to consider outlier. Defaults to 4
% num        - number of outliers a subject can have before being flagged
%
% Outputs
%
% outliers   - binary vector denoting which subjects are outliers
% afq        - if the afq structure is an output then outliers will be
%              appended to afq.metadata.outliers
%
% Example:
%
% [outliers ,afq] = AFQ_outliers(afq, {'dki_MK', 'dki_MD'}, 4, 40)

%% Argument checking
if ~exist('properties', 'var') || isempty(properties)
    properties = {'md'};
end
if ~exist('thresh', 'var') || isempty(thresh)
    thresh = 4;
end
if ~exist('num', 'var') || isempty(num)
    num = 20;
end
%% Computation
outliers = zeros(length(afq.sub_dirs),1);
for ii = 1:length(properties)
    vals = AFQ_get(afq, 'vals', properties{ii});
    m = nanmean(vals);
    s = nanstd(vals);
    zscores = (vals - repmat(m,size(vals,1),1))./repmat(s,size(vals,1),1);
    outliers = outliers + sum(abs(zscores) > thresh ,2);
end

% Find subjects that excede threshold
outliers = outliers >= num;

if nargout == 2
    afq.metadata.outliers = outliers;
end
