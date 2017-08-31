function afq = AFQ_RemoveSubjects(afq, subnums)
% Remove subjects from an AFQ structure
%
% Code to remove subjects from a precomputed afq structure.
%
% Inputs:
%
% afq     - afq structure. See AFQ_run
% subnums - indices of subjects to remove
%
% Outputs:
%
% afq     - afq structure with subjects removed
%
% Example:
%
% load('~/git/lifespan/data/WH_database_full_metadata.mat');
% rmsubs = find(afq.metadata.Age>55);
% afq = AFQ_RemoveSubjects(afq, rmsubs);

valNames = fieldnames(afq.vals);
fgNames  = AFQ_get(afq,'fgnames');
nsubs = length(afq.sub_dirs);
keep = ones(nsubs,1);
keep(subnums) = 0;
keep = find(keep);

% Loop over fiber groups and values to remove for unwanted subjects
for jj = 1:length(fgNames)
    for v = 1:length(valNames)
        % Remove values for unwanted subjects. We are removing the rows from the
        % matrix
        afq.vals.(valNames{v}){jj} = afq.vals.(valNames{v}){jj}(keep,:);
    end
end

% Remove paths to fiber groups
fgFiles = fieldnames(afq.files.fibers);
for jj = 1:length(fgFiles)
    if length(afq.files.fibers.(fgFiles{jj})) == length(afq.sub_group)
        afq.files.fibers.(fgFiles{jj}) = afq.files.fibers.(fgFiles{jj})(keep)
    end
end

% Next get remove Tract Profiles for unwanted subjects
afq.TractProfiles = afq.TractProfiles(keep,:);

% Remove unwanted subject directories
afq.sub_dirs = afq.sub_dirs(keep);
afq.sub_group = afq.sub_group(keep);

% And metadata
mdata = fieldnames(afq.metadata)
for ii = 1:length(mdata)
    afq.metadata.(mdata{ii}) =  afq.metadata.(mdata{ii})(keep);
end