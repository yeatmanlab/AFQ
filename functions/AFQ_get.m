function val = AFQ_get(afq, param, varargin)
% Get value from afq structure
%
% val = AFQ_get(afq, param, varargin)
%
% param list:              - arguments:
%
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
% 'all vals'              - 'valname' , 'group'
% 'use mrtrix'
% 'dt6path'               - [subject number]
% 'tracking parameters'
% 'mrtrixpaths'           - [subject number]
% Written by Jason D. Yeatman August 2012

% remove spaces and upper case
param = mrvParamFormat(param);

switch(param)
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
    case{'vals' 'allvals' 'valsall'}
        if ~exist('varargin' ,'var') || isempty(varargin)
            error('Need to define which value: AFQ_get(afq,''vals'',''fa'')');
        elseif(nargin == 3)
            % First argument is the value name
            valnames = fieldnames(afq.vals);
            v = find(strcmpi(varargin{1},valnames));
            if ~isempty(v)
                valname = valnames{v};
                val = horzcat(afq.vals.(valname){:});
            else
                error('%s values do not exist',varargin{1})
            end
        elseif(nargin == 4)
            % Second argument is the subject group
            switch(varargin{2})
                case{'patient', 1}
                    valnames = fieldnames(afq.patient_data);
                    v = find(strcmpi(varargin{1},valnames));
                    if ~isempty(v)
                        valname = valnames{v};
                        val = horzcat(afq.patient_data(:).(valname));
                    else
                        error('%s values do not exist',varargin{1})
                    end
                case{'control', 0}
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
    otherwise
        error('Uknown afq parameter');
end

return
