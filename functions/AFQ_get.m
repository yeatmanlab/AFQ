function val = AFQ_get(afq, param, varargin)
% Get value from afq structure
%
% val = AFQ_get(afq, param, varargin)
%
% param list:              - arguments:
%
% 'fiber group names'
% 'number of images'
% 'subject group'
% 'patient data'
% 'control data'
% 'norms'
% 'number of patients'
% 'number of controls'
% 'number of subjects'
% 'number of fibergroups'
% 'track whole brain'     - [subject number]
% 'segmentfibers'         - [subject number]
% 'cleanfibers'           - [subject number]
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
% 'use mrtrix'
% 'dt6path'               - [subject number]
% 'tracking parameters'
% 'mrtrixpaths'           - [subject number]
% 'show figures'
% 
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
    case{'fibergroupnames' 'fgnames'}
        val = afq.fgnames;
    case{'numimages', 'numberofimages'}
        val = length(afq.files.images);
    case({'sub_group' 'subgroup' 'subjectgroup'})
        val = afq.sub_group;
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
    case({'numberofsubjects' 'numsubjects' 'numsubs'})
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
    case{'docleaning'}
        % Check if user wants to overwrite cleaned fibers for this subject or
        % Fibers have not yet been cleaned 
        val = logical(afq.overwrite.fibers.clean(varargin{1})) || ...
            isempty(afq.files.fibers.clean{varargin{1}}) || ...
            ~ischar(afq.files.fibers.clean{varargin{1}});
    case{'cleanfibers' 'cleanedfibers' 'cleanfg'}
            val = dtiReadFibers(afq.files.fibers.clean{varargin{1}});
    case{'computevals' 'computeprofiles' 'computetractprofiles' 'compute'}
        % Check if user wants to overwrite values for this subject or
        % No values have been computed yet or
        % Values have not yet been computed for this subject
        val = logical(afq.overwrite.vals(varargin{1})) || ...
            isempty(afq.vals.fa) || ...
            size(afq.vals.fa{1},1) < varargin{1};
    case{'vals' 'allvals' 'valsall' fgnames{:}}
        if ~exist('varargin' ,'var') || isempty(varargin)
            error('Need to define which value: AFQ_get(afq,''vals'',''fa'')');
            % Three arguments means that the user defined the value they
            % want but no group
        elseif(nargin == 3)
            % First argument is the value name
            valnames = fieldnames(afq.vals);
            v = find(strcmpi(varargin{1},valnames));
            % Check if the user defined a specific fiber group
            fgNum = find(strcmpi(param,fgnames));
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
                        val = horzcat(afq.patient_data(:).(valname));
                    else
                        error('%s values do not exist',varargin{1})
                    end
                case{'control','controls', 0}
                    valnames = fieldnames(afq.control_data);
                    v = find(strcmpi(varargin{1},valnames));
                    if ~isempty(v)
                        valname = valnames{strcmpi(varargin{1},valnames)};
                        val = horzcat(afq.control_data(:).(valname));
                    else
                        error('%s values do not exist',varargin{1})
                    end
                otherwise
                    error('Do you want vals for patients or controls?')
            end
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
        val = afq.files.dt6{varargin{1}};
    case{'trackingparameters'}
        val = afq.params.track;
    case{'mrtrixpath' 'mrtrixpaths'}
        val.csd = afq.files.mrtrix.csd{varargin{1}};
        val.wm  = afq.files.mrtrix.wm{varargin{1}};
    case{'showfigures' 'showfigs'}
        val = logical(afq.params.showfigs);
    case{'fiberweighting'}
        val = afq.params.fiberWeighting;
    otherwise
        error('Uknown afq parameter');
end

return
