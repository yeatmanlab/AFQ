function afq = AFQ_SubjectAvgMetadata(afq, subIDs)
% Create new metadata variables that contain ave values for each subject
%
% afq = AFQ_SubjectAvgMetadata(afq, subIDs)
%
% This function is useful for longitudinal data where you might want to
% create average variables within subjects
%
% Example:
%
% afq = afq = AFQ_SubjectAvgMetadata(afq, afq.sub_names)

if ~exist('subIDs', 'var') || isempty(subIDs) && exist(afq.sub_names)
    subIDs = afq.sub_names;
end

varnames = fieldnames(afq.metadata);
us = unique(subIDs);
for ii = 1:numel(varnames)
    if isnumeric(afq.metadata.(varnames{ii}))
        ssvar = nan(size(afq.metadata.(varnames{ii})));
        ssdmvar = nan(size(afq.metadata.(varnames{ii})));
        for ss = 1:numel(us)
            idx = strcmp(subIDs,us{ss});
            tmp = afq.metadata.(varnames{ii})(idx);
            ssvar(idx) = repmat(mean(tmp),size(tmp));
            ssdmvar(idx) = tmp-mean(tmp);
        end
        afq.metadata.(['sm_' varnames{ii}]) = ssvar; clear ssvar;
        afq.metadata.(['sdm_' varnames{ii}]) = ssdmvar; clear ssdmvar;
    end
end