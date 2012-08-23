function afq = AFQ_set(afq,param,varargin)
% Set properties of afq object
%
% afq = AFQ_set(afq, param, varargin)
%
% Parameter list and associated values:
%
%  param     - varargin
%
% 'images'   - 1xN cell array (N is number of subjects) of paths to images
%              to perform AFQ computations on
% 'vals'     - 'subnum', subnum, 'valnam', vals
% 'norms'    -
% 'sub_group'- sub_group

% remove spaces and upper case
param = mrvParamFormat(param);

switch(param)
    case 'images'
        afq.files.images(end+1).path = varargin{1};
        if length(varargin) > 1
            afq.files.images(end+1).name = varargin;
        else
            [path, name, ext] = fileparts(val);
            afq.files.images(end+1).name = name;
        end
        
    case 'vals'
        if length(varargin) >= 4
            subnum = varargin{2};
            for ii = 3:2:length(varargin)
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
        
    case 'norms'
        valnames = fieldnames(afq.vals);
        for ii = 1:length(valnames)
            for jj = 1:length(afq.fgnames)
                afq.norms.(['mean' upper(valnames{ii})])(:,jj)  = nanmean(afq.vals.(valnames{ii}){jj}(afq.sub_group==0,:));
                afq.norms.(['sd' upper(valnames{ii})])(:,jj)  = nanstd(afq.vals.(valnames{ii}){jj}(afq.sub_group==0,:));
                afq.patient_data(jj).([upper(valnames{ii})]) = afq.vals.(valnames{ii}){jj}(afq.sub_group==1,:);
                afq.control_data(jj).([upper(valnames{ii})]) = afq.vals.(valnames{ii}){jj}(afq.sub_group==0,:);
            end
        end
        
    case 'sub_group'
        afq.sub_group = varargin{1};
        
    case 'currentsubject'
        afq.currentsub = varargin{1};
end



