function [outliers, afq] = AFQ_outliers(afq, properties, thresh, num, mexclude, dvols)
% Detect outliers who have tracts with unexpected values
%
% outliers = AFQ_outliers(afq, properties, thresh, num, mexclude)
%
% Inputs
%cd 
% afq        - afq structure. See AFQ_run
% properties - properties upon which to detect outliers. Cell array
% thresh     - distance (z-score) to consider outlier. Defaults to 4
% num        - number of outliers a subject can have before being flagged
% mexclude   - binary. exlude datasets with excessive motion (1) or not
% dvols      - binary. exlude datasets with excessive # of dropped volumes
%
% Outputs
%
% outliers   - binary vector denoting which subjects are outliers
% afq        - if the afq structure is an output then outliers will be
%              appended to afq.metadata.outliers
%
% Example:
%
% [outliers ,afq] = AFQ_outliers(afq, {'dki_MK', 'dki_MD'}, 4, 40, 1)

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
if ~exist('mexclude', 'var') || isempty(mexclude)
    mexclude = 0; % don't exclude based on motion by default
end
if ~exist('dvols', 'var') || isempty(mexclude)
    dvols = 0; % don't exclude based on # of dropped volumes
end

%% Computation
% First exclude subjects/sessions with values outside a
% reasonable range (defined by 'thresh'):
outliers = zeros(length(afq.sub_dirs),1);
for ii = 1:length(properties)
    vals = AFQ_get(afq, 'vals', properties{ii});
    m = nanmean(vals);
    s = nanstd(vals);
    zscores = (vals - repmat(m,size(vals,1),1))./repmat(s,size(vals,1),1);
    outliers = outliers + sum(abs(zscores) > thresh, 2);
end
% Find subjects that excede threshold
outliers = outliers >= num;

% Next, optionally define outliers based on motion:
if mexclude == 1
    motion = AFQ_get(afq,'metadata','motion');
    motion = motion > 0.7;
    % exlude subjects based on motion
    outliers(motion(1:length(outliers)) == 1) = 1;
end

% Optionally exclude subjects with > 10% dropped volumes per scan:
if dvols == 1
    totalvols = 111; % total number of volumes in our diffusion seq
    dropped = discardedVolumes(afq.sub_names, afq.metadata.session);
    outliers((dropped./totalvols) > 0.1) = 1;
end

if nargout == 2
    afq.metadata.outliers = outliers;
end


