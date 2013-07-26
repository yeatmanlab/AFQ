function afq = AFQ_CombineSgeRuns(afq, afq_sge, subnum)
% Add the values from an individual subject's sge run to an afq structure
%
% afq = AFQ_CombineSgeRuns(afq, afq_sge, subnum)

valNames = fieldnames(afq.vals);
fgNames  = AFQ_get(afq,'fgnames');

% Loop over fiber groups and values to assign the values for this subjects
% sge run to the main afq structure
for jj = 1:length(fgNames)
   for v = 1:length(valNames)
       % Set the values for this subject in the afq struct to be the ones
       % from the afq_sge struct for this valname and fiber group
       afq.vals.(valNames{v}){jj}(subnum,:) = afq_sge.vals.(valNames{v}){jj}(subnum,:); 
   end
end

% Next get that subject's TractProfiles and add them to the main afq struct
afq.TractProfiles(subnum,:) = afq_sge.TractProfiles(subnum,:);