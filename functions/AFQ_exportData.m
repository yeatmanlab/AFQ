function AFQ_exportData(data, filename, property)
% Export AFQ results as text files that can be read into excel or spss.
% 
% AFQ_exportData(data, filename, property)
%
% Inputs:
% data     = data output from AFQ
% filename = name of the file to write including the path
% property = a string denoting which diffusion property to plot. 
%             options are 'fa' 'rd' 'ad' 'md'
%
% Example:
%
% [patient_data control_data] = AFQ_run(sub_dirs, sub_group); % See AFQ_run
% AFQ_exportData(control_data, '/home/jyeatman/results', 'fa');
%
% (c) Jason D. Yeatman, Vista Team, 4/3/2012
%
if ~exist('property','var') || isempty(property)
    property = 'fa';
else
    property = lower(property); %make sure property is lower case
end
% These are the names of the fiber groups
fgNames={'Left Thalmic Radiation','Right Thalmic Radiation','Left Corticospinal','Right Corticospinal', 'Left Cingulum Cingulate', 'Right Cingulum Cingulate'...
    'Left Cingulum Hippocampus','Right Cingulum Hippocampus', 'Callosum Forceps Major', 'Callosum Forceps Minor'...
    'Left IFOF','Right IFOF','Left ILF','Right ILF','Left SLF','Right SLF','Left Uncinate','Right Uncinate','Left Arcuate','Right Arcuate'};
% write out data for each fiber group as a separate file
switch(property)
    case 'fa'
        for jj = 1:length(data)
            csvwrite([filename '_' fgNames{jj} '_' property], data(jj).FA);
        end
    case 'rd'
        for jj = 1:length(data)
            csvwrite([filename '_' fgNames{jj} '_' property], data(jj).RD);
        end
    case 'ad'
        for jj = 1:length(data)
            csvwrite([filename '_' fgNames{jj} '_' property], data(jj).AD);
        end
    case 'md'
        for jj = 1:length(data)
            csvwrite([filename '_' fgNames{jj} '_' property], data(jj).MD);
        end
end
