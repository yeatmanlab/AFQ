function afq = AFQ_set(afq,param,varargin)
% Set properties of afq object
%
% afq = AFQ_set(afq, param, [vargrgin])
%
% The afq structure stores all the AFQ computations.  AFQ_set is the main
% function to add values, images, computations etc. to this structure.  As
% inputs AFQ_set takes in an afq structure, a parameter defined as a
% 'string' which is to be set, and a number of aditional arguments and
% values that are specific to that parameter.  The arguments for each
% parameter are described below.  AFQ_set always returns the modified afq
% structure.  The structure can be queried and values can be accessed with
% AFQ_get.  The structure is created with AFQ_Create.
%
% See also: AFQ_get AFQ_Create
%
% Parameter list and associated arguments:
%
%  param           - varargin
%
% 'images'         - Add images to compute tract profiles on.
%                    varargin = 1xN cell array (N is number of subjects) of
%                    paths to nifti images.
% 'vals'           - Add Tract Profile values to afq structure.
%                    varargin = 'subnum', subnum, 'valname', vals
% 'norms'          - Compute and assign control group norms to afq.norms
%                    varargin = no argument needed
% 'sub_group'      - Define subject group (patient=1 control=0) to afq
%                    structure.
%                    varargin = [1 1 1 0 0 0.....,]
% 'current subject'- Define current subject for afq computations
%                    varargin = [subject number]
% 'overwritefibers'- Recompute fibers for a subject. If the second argument
%                    is blank then it will recompute for all subjects
%                    varargin = [subject number]
% 'overwritevals'  - Recompute Tract Profile values for a subject
%                    varargin = [subject number]
% 'wholebrain fg path'- Set the path to the wholebrain fiber group for a
%                       subject
%                       varargin = 'subnum', [subject number]
% 'segmented fg path' - Set the path to the segmented fiber group for a
%                       subject
%                       varargin = 'subnum', [subject number]
% 'clean fg path'     - Set the path to the cleaned fiber group for a
%                       subject
%                       varargin = 'subnum', [subject number]
% 'tractprofiles'  - Add tract profiles to the afq structure for a subject.
%                    This should be a structure with each of the AFQ fiber
%                    groups in it.  See AFQ_CreateTractProfile
%                    varargin = 'subnum', [subject number], TractProfiles
%
% Example:
%
% afq = AFQ_set(afq, 'segemented fg path','subnum',10,'/data/sub10/fibers/'Morigroups.mat')
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
                % The first value will be in varargin{3}
                val1 = 3;
            else
                % Otherwise add to the row of the current subject that is
                % being processed
                subnum = afq.currentsub;
                % The first value will be in varargin{1}
                val1 = 1;
            end
            % Loop over the values that were input
            for ii = val1:2:length(varargin)
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
        
    case 'sub_group' % Define subject groups in afq.sub_group
        afq.sub_group = varargin{1};
        
    case 'currentsubject' % Define the current subject being analyzed
        afq.currentsub = varargin{1};
        
    case 'overwritefibers' % Recompute fiber groups for subject # varargin
        % If a subject number is defined only overwrite for that subject
        if nargin == 3
            afq.overwrite.fibers.wholebrain(varargin{1}) = 1;
            afq.overwrite.fibers.segmented(varargin{1})  = 1;
            afq.overwrite.fibers.clean(varargin{1})      = 1;
        else
            % Otherwise overwrite for all subjects
            afq.overwrite.fibers.wholebrain(:)= 1;
            afq.overwrite.fibers.segmented(:)= 1;
            afq.overwrite.fibers.clean(:)= 1;
        end
        
    case 'overwritevals' % Recompute values for subject # varargin
        afq.overwrite.vals(varargin{1}) = 1;
        
    case {'wholebrainfgpath' 'wholebrainpath' 'wholebrainfg' 'wholebrainfibergroup'}
        % Add the values to the correct row of the data matrix
        if strcmp('subnum',varargin{1})
            % If the subject number was defined use that row
            subnum = varargin{2};
            % The first value will be in varargin{3}
            val1 = 3;
        else
            % Otherwise add to the row of the current subject that is
            % being processed
            subnum = afq.currentsub;
            % The first value will be in varargin{1}
            val1 = 1;
        end
        afq.files.fibers.wholebrain{subnum} = varargin{val1};
    case {'segmentedfgpath' 'segmentedpath' 'segmentedfg' 'segmentedfibergroup'}
        % Add the values to the correct row of the data matrix
        if strcmp('subnum',varargin{1})
            % If the subject number was defined use that row
            subnum = varargin{2};
            % The first value will be in varargin{3}
            val1 = 3;
        else
            % Otherwise add to the row of the current subject that is
            % being processed
            subnum = afq.currentsub;
            % The first value will be in varargin{1}
            val1 = 1;
        end
        afq.files.fibers.segmented{subnum} = varargin{val1};
    case {'cleanfgpath' 'cleanpath' 'cleanfg' 'cleanfibergroup'}
        % Add the values to the correct row of the data matrix
        if strcmp('subnum',varargin{1})
            % If the subject number was defined use that row
            subnum = varargin{2};
            % The first value will be in varargin{3}
            val1 = 3;
        else
            % Otherwise add to the row of the current subject that is
            % being processed
            subnum = afq.currentsub;
            % The first value will be in varargin{1}
            val1 = 1;
        end
        afq.files.fibers.clean{subnum} = varargin{val1};
    case{'tractprofiles' 'tractprofile'}
        % Add the values to the correct row of the data matrix
        if strcmp('subnum',varargin{1})
            % If the subject number was defined use that row
            subnum = varargin{2};
            % The first value will be in varargin{3}
            val1 = 3;
        else
            % Otherwise add to the row of the current subject that is
            % being processed
            subnum = afq.currentsub;
            % The first value will be in varargin{1}
            val1 = 1;
        end
        % Loop over the number of tract profiles and assign them to the afq
        % structure
        for jj = 1:length(varargin{val1})
            afq.TractProfiles(subnum,jj) = varargin{val1}(jj);
        end
    otherwise
        error('Uknown afq parameter');
end



