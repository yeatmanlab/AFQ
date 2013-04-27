function val = AFQ_get(afq, param, varargin)
% Get value from afq structure
%
% val = AFQ_get(afq, param, varargin)
%
% param list:              - arguments:
%
% 'sub_dirs'
% 'fiber group names'
% 'number of images'
% 'subject group'
% 'patient data'
% 'control data'
% 'norms'
% 'numfg'
% 'number of patients'
% 'number of controls'
% 'number of subjects'
% 'number of fibergroups'
% 'track whole brain'     - [subject number]
% 'wholebrain fg'         - [subject number]
% 'segmentfibers'         - [subject number]
% 'cleanfibers'           - [subject number]
% 'computenorms' 
%
% % To get the path to a fiber group
% 'fgname path'           - [subject number]
%
% 'computevals'           - [subject number]
%
% To get all the values for a particular measurement (eg. 'fa') for all the
% subjects or a group of subjects (eg. 'patients'). If no input is passed
% in for the group then it will return the values for all the subjects.
% 'all vals'              - 'valname' , ['group']
%
% To get the values for a particular measurement (eg. 'fa') for a particlar
% tract for all the subjects or a group of subjects. For the different
% fiber group names see AFQ_get(afq, 'fgnames')
% 'tractname'             - 'valname', 'group'
% 'left arcuate'          - 'fa'     , 'controls'
%
% To get values that are saved within the tract profiles rather than in the
% vals field:
% 'TractProfile vals'      - 'tract name', 'valname'
% 'TractProfile vals'      - 'Left Arcuate', 'fa'
%
% 'use mrtrix'
% 'dt6path'               - [subject number]
% 'tracking parameters'
% 'mrtrixpaths'           - [subject number]
% 'show figures'
% 'roi1 path'             - [fg number], [subject number]
% 'roi2 path'             - [fg number], [subject number]
% 'current subject'
% 'seg fg name'           - [subject number]
% 'clean fg name'         - [subject number]
% 'segfilename'            -
% To get any of the parameters save in the afq structure (see AFQ_Create),
% enter the name of the parameter. Some have not been implimented yet, but
% will be soon.
% AFQ_get(afq,'fiberWeighting')
%
% Written by Jason D. Yeatman August 2012

% remove spaces and upper case
param = mrvParamFormat(param);
% Get the names of all the fiber groups without spaces or capitals
fgnames = cell(1,length(afq.fgnames));
for jj = 1:length(fgnames)
    fgnames{jj} = mrvParamFormat(afq.fgnames{jj});
end

%% Get the requested parameter
switch(param)
    case{'sub_dirs' 'subs' 'subjectdirectories'}
        val = afq.sub_dirs;
        if isempty(val)
            error('no subject directories');
        end
    case{'fibergroupnames' 'fgnames'}
        val = afq.fgnames;
    case{'numimages', 'numberofimages'}
        val = length(afq.files.images);
    case({'sub_group' 'subgroup' 'subjectgroup'})
        val = afq.sub_group;
        if isempty(val)
            error('no subject group');
        end
    case({'patient_data' 'patientdata'})
        val = afq.patient_data;
    case({'control_data' 'controldata'})
        val = afq.control_data;
    case('norms')
        val = afq.norms;
    case({'numberofpatients' 'numpatients'})
        val = sum(afq.sub_group);
    case({'numberofcontrols' 'numcontrols'});
        val = sum(afq.sub_group == 0);
    case({'numberofsubjects' 'numsubjects' 'numsubs' 'nsubs'})
        val = length(afq.sub_group);
    case({'numberoffibergroups' 'numfg' 'numfibergroups' 'nfg' 'numberfibergroups'});
        val = length(afq.fgnames);
    case{'dotracking' 'trackfibers' 'track'}
        % Check if user wants to overwrite wholebrain tractography or
        % Wholebrain fiber group has not been tracked
        val = logical(afq.overwrite.fibers.wholebrain(varargin{1})) || ...
            isempty(afq.files.fibers.wholebrain{varargin{1}}) || ...
            ~ischar(afq.files.fibers.wholebrain{varargin{1}});
    case{'wholebrainfibergroup' 'wholebrain' 'wholebrainfg'}
        val = dtiReadFibers(afq.files.fibers.wholebrain{varargin{1}});
    case{'dosegmentation'}
        % Check if user wants to overwrite segmented fibers or
        % Fibers have not yet been segmented
        val = logical(afq.overwrite.fibers.segmented(varargin{1})) || ...
            isempty(afq.files.fibers.segmented{varargin{1}}) || ...
            ~ischar(afq.files.fibers.segmented{varargin{1}});
    case{'segmentedfibers' 'morigroups'}
        val = dtiReadFibers(afq.files.fibers.segmented{varargin{1}});
    case horzcat(strcat(fgnames,'path'))
        % find the fiber group number
        n = find(strcmpi(param(1:end-4),fgnames));
        if n <= 20
            % If it is one of the mori groups than get the path to that
            % fiber group (preferably cleaned)
            try
                val = afq.files.fibers.clean{varargin{1}};
            catch
                val = afq.files.fibers.segmented{varargin{1}};
            end
        else
            % get the name (because we formated the parameter)
            name = afq.fgnames{n};
            try
                val = afq.files.fibers.([name 'clean']){varargin{1}};
            catch
                val = afq.files.fibers.(name){varargin{1}};
            end
        end
    case horzcat(strcat(fgnames,'fg'))
        % find the fiber group number (get rid of the 'fg' at the end
        n = find(strcmpi(param(1:end-2),fgnames));
        if n <= 20
            % If it is one of the mori groups than get the path to that
            % fiber group (preferably cleaned)
            try
                path = afq.files.fibers.clean{varargin{1}};
            catch
                path = afq.files.fibers.segmented{varargin{1}};
            end
            % Load the fiber group
            fg = dtiReadFibers(path);
            % Pull out the desired fiber group number
            val = fg(n);
        else
            % get the name (because we formated the parameter)
            name = afq.fgnames{n};
            try
                path = afq.files.fibers.([name 'clean']){varargin{1}};
            catch
                path = afq.files.fibers.(name){varargin{1}};
            end
            % Load the fiber group
            val = dtiLoadFiberGroup(path);
        end
    case{'docleaning'}
        % Check if user wants to overwrite cleaned fibers for this subject or
        % Fibers have not yet been cleaned
        val = logical(afq.overwrite.fibers.clean(varargin{1})) || ...
            isempty(afq.files.fibers.clean{varargin{1}}) || ...
            ~ischar(afq.files.fibers.clean{varargin{1}});
    case{'cleanfibers' 'cleanedfibers' 'cleanfg'}
        val = dtiReadFibers(afq.files.fibers.clean{varargin{1}});
    case{'segname' 'segmentedfibersname' 'segfgname'}
        [~, val] = fileparts(afq.files.fibers.segmented{varargin{1}});
    case{'cleanfgname'}
        [~, val] = fileparts(afq.files.fibers.clean{varargin{1}});
    case{'computevals' 'computeprofiles' 'computetractprofiles' 'compute'}
        % Check if user wants to overwrite values for this subject or
        % No values have been computed yet or
        % Values have not yet been computed for this subject
        val = logical(afq.overwrite.vals(varargin{1})) || ...
            isempty(afq.vals.fa) || ...
            size(afq.vals.fa{1},1) < varargin{1};
    case 'computenorms'
        % Check if norms should be computed
        if isfield(afq.params,'computenorms') && ~isempty(afq.params.computenorms)
            val = logical(afq.params.computenorms);
        else
            val = true;
        end
    case{'vals' 'allvals' 'valsall' fgnames{:}}
        % Check if the user defined a specific fiber group and if so get
        % that fiber group number
        fgNum = find(strcmpi(param,fgnames));
        if ~exist('varargin' ,'var') || isempty(varargin)
            error('Need to define which value: AFQ_get(afq,''vals'',''fa'')');
            % Three arguments means that the user defined the value they
            % want but no group
        elseif(nargin == 3)
            % First argument is the value name
            valnames = fieldnames(afq.vals);
            v = find(strcmpi(varargin{1},valnames));
            if ~isempty(v)
                valname = valnames{v};
                % If a specific fiber group was defined get vals for that
                % group
                if ~isempty(fgNum) && length(afq.vals.(valname)) >= fgNum
                    val = afq.vals.(valname){fgNum};
                elseif ~isempty(fgNum)
                    error('%s values do not exist',fgnames{fgNum});
                else
                    % Otherwise get vals for all groups
                    val = horzcat(afq.vals.(valname){:});
                end
            else
                error('%s values do not exist',varargin{1})
            end
            % Four arguments means the user also defined the group
        elseif(nargin == 4)
            % Second argument is the subject group
            switch(varargin{2})
                case{'patient','patients', 1}
                    valnames = fieldnames(afq.patient_data);
                    v = find(strcmpi(varargin{1},valnames));
                    if ~isempty(v)
                        valname = valnames{v};
                        % If a specific fiber group was defined get vals for that
                        % group
                        if ~isempty(fgNum) && length(afq.patient_data) >= fgNum
                            val = afq.patient_data(fgNum).(valname);
                        elseif ~isempty(fgNum)
                            error('%s values do not exist',fgnames{fgNum});
                        else
                            % Otherwise get vals for all groups
                            val = horzcat(afq.patient_data(:).(valname));
                        end
                    else
                        error('%s values do not exist',varargin{1})
                    end
                case{'control','controls', 0}
                    valnames = fieldnames(afq.control_data);
                    v = find(strcmpi(varargin{1},valnames));
                    if ~isempty(v)
                        valname = valnames{strcmpi(varargin{1},valnames)};
                        % If a specific fiber group was defined get vals for that
                        % group
                        if ~isempty(fgNum) && length(afq.control_data) >= fgNum
                            val = afq.control_data(fgNum).(valname);
                        elseif ~isempty(fgNum)
                            error('%s values do not exist',fgnames{fgNum});
                        else
                            % Otherwise get vals for all groups
                            val = horzcat(afq.control_data(:).(valname));
                        end
                    else
                        error('%s values do not exist',varargin{1})
                    end
                otherwise
                    error('Do you want vals for patients or controls?')
            end
        end
    case {'tractprofilevals' 'valstractprofile'}
        % Get the number of the tract. This is the first input argumet
        tnum = find(strcmp(varargin{1},AFQ_get(afq,'fgnames')));
        if isempty(tnum)
            fprintf('\n please provide valid tract name\n')
            fprintf('correct call AFQ_get(afq,''tractprofilevals'',''Left Arcuate'',''fa'')\n')
        end
        % Loop over each subject's tract profile and get the requested
        % values
        valname = varargin{2};
        for ii = 1:AFQ_get(afq, 'numsubs')
           val(ii,:) = afq.TractProfiles(ii,tnum).vals.(valname);
        end
    case{'usemrtrix'}
        if afq.software.mrtrix == 1 && ...
                strcmp('mrtrix',afq.params.track.algorithm)...
                && afq.params.computeCSD == 1;
            val = true;
        else
            val = false;
        end
    case{'dt6path'}
        % If a subject number was input then get the dt6path for that
        % subject
        if nargin == 3
            val = afq.files.dt6{varargin{1}};
            % Otherwise return a cell array of paths for all subjects
        else
            val  = afq.files.dt6;
        end
    case{'trackingparameters'}
        val = afq.params.track;
    case{'mrtrixpath' 'mrtrixpaths'}
        val.csd = afq.files.mrtrix.csd{varargin{1}};
        val.wm  = afq.files.mrtrix.wm{varargin{1}};
    case{'showfigures' 'showfigs'}
        val = logical(afq.params.showfigs);
    case{'fiberweighting'}
        val = afq.params.fiberWeighting;
    case{'clip2rois'}
        val = afq.params.clip2rois;
    case {'roi1' 'roi1path'}
        if nargin~=4
            error('correct call: AFQ_get(afq, ''roi1'',[roi number], [subject number])')
        else
            % varargin{1} is the roi number varargin{2} is the subject number
            val = fullfile(afq.sub_dirs{varargin{2}},'ROIs',afq.roi1names{varargin{1}});
        end
    case {'roi2' 'roi2path'}
        if nargin~=4
            error('correct call: AFQ_get(afq, ''roi1'',[roi number], [subject number])')
        else
            % varargin{1} is the roi number varargin{2} is the subject number
            val = fullfile(afq.sub_dirs{varargin{2}},'ROIs',afq.roi2names{varargin{1}});
        end
    case {'currentsubject' 'cursub'}
        val = afq.currentsub;
    case {'outdir' 'outputdirectory'}
        val = afq.params.outdir;
    case {'outname' 'outputname'}
        if isfield(afq.params,'outname')
            val = afq.params.outname;
        else
            val = [];
        end
    case {'runsubs' 'runsubjects'}
        if ~isfield(afq, 'runsubs') || isempty(afq.runsubs) || ischar(afq.runsubs)
            val = 1:AFQ_get(afq,'numberofsubjects');
        else
            val = afq.runsubs;
        end
        % transpose if it's a column vector
        if size(val,1) > size(val,2)
            val = val';
        end
    case 'segfilename'
        if isfield(afq.files.fibers,'segName')
            val = afq.files.fibers.segName;
        elseif AFQ_get(afq,'clip2rois') == 1
            val = 'MoriGroups.mat';
        elseif AFQ_get(afq,'clip2rois') == 0
            val = 'MoriGroups_Cortex.mat';
        end
    otherwise
        error('Uknown afq parameter');
end

return
