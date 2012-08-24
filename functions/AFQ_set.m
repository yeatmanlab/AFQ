function afq = AFQ_set(afq,param,varargin)
% Set properties of afq object
%
% afq = AFQ_set(afq, param, varargin)
%
% Parameter list and associated values:
%
%  param     - varargin
%
% 'images'         - Add images to compute tract profiles on. 
%                    varargin = 1xN cell array (N is number of subjects) of
%                    paths to nifti images.
% 'vals'           - Add Tract Profile values to afq structure. 
%                    varargin = ('subnum', subnum, 'valnam', vals 'norms')         -
% 'sub_group'      - Define subject group (patient=1 control=0) to afq 
%                    structure.  
%                    varargin = [1 1 1 0 0 0.....,] 
% 'currentsubject' - Define current subject for afq computations
%                    varargin = subjectnumber (eg. 10)
% 'norms'          - Compute and assign control group norms to afq.norms
%                    varargin = no argument needed
%
% Written by Jason D. Yeatman August 2012

% remove spaces and upper case
param = mrvParamFormat(param);

switch(param)
    
    case 'images' % Add images to afq.files.images
        afq.files.images(end+1).path = varargin{1};
        if length(varargin) > 1
            afq.files.images(end).name = varargin{2};
        else
            s = strfind(varargin{1}{1},'/');
            p = strfind(varargin{1}{1},'.');
            p = p(p > s(end));
            name = varargin{1}{1}(s(end)+1:p(1)-1);
            afq.files.images(end).name = name;
        end
        
    case 'vals' % Add values to the afq.values field
        if length(varargin) >= 4
            % Add the values to the correct row of the data matrix
            if strcmp('subnum',varargin{1})
                % If the subject number was defined use that row
                subnum = varargin{2};
            else
                % Otherwise add to the row of the current subject that is
                % being processed
                subnum = afq.currentsub;
            end
            % Loop over the values that wer input
            for ii = 3:2:length(varargin)
                % Loop over the fiber tracts
                for jj = 1:20
                    % Take the stats that were calculated in
                    % AFQ_ComputeTractProperties and add them to a sructure
                    % for the full sample of subjects.  Each fiber group
                    % has its own cell. Within each cell there is a row for
                    % each subject with numberofnodes columns
                    afq.vals.(varargin{ii}){jj}(subnum,:) = varargin{ii + 1}(:, jj);
                end
            end
        end
        
    case 'norms' % Compute and assign control group norms to afq.norms
        valnames = fieldnames(afq.vals);
        for ii = 1:length(valnames)
            for jj = 1:length(afq.fgnames)
                afq.norms.(['mean' upper(valnames{ii})])(:,jj)  = nanmean(afq.vals.(valnames{ii}){jj}(afq.sub_group==0,:));
                afq.norms.(['sd' upper(valnames{ii})])(:,jj)  = nanstd(afq.vals.(valnames{ii}){jj}(afq.sub_group==0,:));
                afq.patient_data(jj).([upper(valnames{ii})]) = afq.vals.(valnames{ii}){jj}(afq.sub_group==1,:);
                afq.control_data(jj).([upper(valnames{ii})]) = afq.vals.(valnames{ii}){jj}(afq.sub_group==0,:);
            end
        end
        
    case 'sub_group' % Add subject grouping to afq.sub_group
        afq.sub_group = varargin{1};
        
    case 'currentsubject' % Define the current subject being analyzed
        afq.currentsub = varargin{1};
end



