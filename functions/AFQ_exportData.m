function AFQ_exportData(data, filename, property, sub_ids)
% Export AFQ results as text files that can be read into excel or spss.
% 
% AFQ_exportData(data, filename, property)
%
% Inputs:
% data     = data output from AFQ
% filename = name of the file to write including the path
% property = a string denoting which diffusion property to plot. 
%             options are 'fa' 'rd' 'ad' 'md'
% sub_ids  = a cell array of subject names or a vector of subject numbers
%
% Example:
%
% [patient_data control_data] = AFQ_run(sub_dirs, sub_group); % See AFQ_run
% % Pull out the subject ids for the controls
% sub_ids = sub_dirs(logical(sub_group));
% % Export the data as a .csv to the AFQ data directory
% [AFQbase AFQdata] = AFQ_directories;
% AFQ_exportData(control_data, [AFQdata '/results'], 'fa', sub_ids);
%
% (c) Jason D. Yeatman, Vista Team, 4/3/2012
%

%% Argument checking

if ~exist('property','var') || isempty(property)
    property = 'FA';
    
    % Check that the property was defined in the right case
elseif isfield(data,lower(property))
    property = lower(property); 
elseif isfield(data,upper(property))
    property = upper(property);
end

if exist('sub_ids','var') && ~isempty(sub_ids)
    % Make sub_ids a column
    if size(sub_ids,2) > size(sub_ids,1)
        sub_ids = sub_ids';
    end
    % Make sure there are as many sub_ids as there are rows in data
    if size(sub_ids,1) ~= size(data(1).(property),1)
        error('There must be a subject ID for each row of data')
    end
    
end

%% Writing the .csv files

if isafq(data)
    fgNames = AFQ_get(data,'fgnames');
end
% These are the names of the fiber groups
fgNames={'Left_Thalmic_Radiation','Right_Thalmic_Radiation','Left_Corticospinal','Right_Corticospinal', 'Left_Cingulum_Cingulate', 'Right_Cingulum_Cingulate'...
    'Left_Cingulum_Hippocampus','Right_Cingulum_Hippocampus', 'Callosum_Forceps_Major', 'Callosum_Forceps_Minor'...
    'Left_IFOF','Right_IFOF','Left_ILF','Right_ILF','Left_SLF','Right_SLF','Left_Uncinate','Right_Uncinate','Left_Arcuate','Right_Arcuate'};
% write out data for each fiber group as a separate file
for jj = 1:length(fgNames)
    % Convert the data for the desired property into a cell array
    if isafq(data)
        m = num2cell(data.vals.(property){jj});
    else
        m = num2cell(data(jj).(property));
    end
    % Add in the subject ids as the first column if they were supplied
    if exist('sub_ids','var') && ~isempty(sub_ids)
        m = horzcat(sub_ids,m);
    end
    % Write out the data file
    cell2csv([filename '_' fgNames{jj} '_' property '.csv'], m);
end

return

%% Origional less efficient code
% switch(property)
%     case 'fa'
%         for jj = 1:length(data)
%             csvwrite([filename '_' fgNames{jj} '_' property], data(jj).FA);
%         end
%     case 'rd'
%         for jj = 1:length(data)
%             csvwrite([filename '_' fgNames{jj} '_' property], data(jj).RD);
%         end
%     case 'ad'
%         for jj = 1:length(data)
%             csvwrite([filename '_' fgNames{jj} '_' property], data(jj).AD);
%         end
%     case 'md'
%         for jj = 1:length(data)
%             csvwrite([filename '_' fgNames{jj} '_' property], data(jj).MD);
%         end
% end
