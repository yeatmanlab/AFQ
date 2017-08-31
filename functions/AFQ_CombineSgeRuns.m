function afq = AFQ_CombineSgeRuns(afq, afq_sge, subnum, subnum_sge)
% Add the values from an individual subject's sge run to an afq structure
%
% afq = AFQ_CombineSgeRuns(afq, afq_sge, subnum, [subnum_sge])
%
% This will add the subject if afq_sge to the afq structure. The default is
% to add subnum to the same place. However, if subnum_sge is defined and
% different from subnum, it will take subnum_sge from afq_sge and add it to
% position subnum with afq. Note that it will only add values from afq_sge
% if there is a corrsponding value in afq.
%
% Input
%
% afq       - afq structure that will contain the data from the new subject
% afq_sge   - afq structure from sge run or from another run of afq
% subnum    - subject number
% subnum_sge- subject number in afq_sge if different from subnum
%
% Output
%
% afq - combined structure

if ~exist('subnum_sge' ,'var') || isempty(subnum_sge)
    subnum_sge = sugnum;
end

% Get the names of the values in the afq_sge structure
valNames = fieldnames(afq_sge.vals);
% Remove values that do not exist in afq structure
indx = ismember(valNames,fieldnames(afq.vals));
valNames = valNames(indx);

fgNames  = AFQ_get(afq_sge,'fgnames');

% Loop over fiber groups and values to assign the values for this subjects
% sge run to the main afq structure
for jj = 1:length(fgNames)
    afq.sub_dirs{subnum} = afq_sge.sub_dirs{subnum_sge};
    for v = 1:length(valNames)
        % Set the values for this subject in the afq struct to be the ones
        % from the afq_sge struct for this valname and fiber group
        afq.vals.(valNames{v}){jj}(subnum,:) = afq_sge.vals.(valNames{v}){jj}(subnum_sge,:);
    end
end
% Check to make sure all the paths to fiber groups are set properly
fgFiles = fieldnames(afq.files.fibers);
for jj = 1:length(fgFiles)
    % Check if this field contains a cell array of paths for each subject
    if iscell(afq_sge.files.fibers.(fgFiles{jj}))
        afq.files.fibers.(fgFiles{jj}){subnum} = afq_sge.files.fibers.(fgFiles{jj}){subnum_sge};
    end
end
% Next get that subject's TractProfiles and add them to the main afq struct
sz = size(afq_sge.TractProfiles);
afq.TractProfiles(subnum,1:sz(2)) = afq_sge.TractProfiles(subnum_sge,1:sz(2));

return