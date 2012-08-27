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
%
%
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
        val = length(afq.fgNames);
    case{'trackwholebrain' 'trackfibers' 'track'}
        val = logical(afq.overwrite.fibers.wholebrain(varargin{1})) || ...
            isempty(afq.files.fibers.wholebrain{varargin{1}}) || ...
            ~ischar(afq.files.fibers.wholebrain{varargin{1}});
    case{'segmentfibers' 'segmentwholebrain'}
        val = logical(afq.overwrite.fibers.segmented(varargin{1})) || ...
            isempty(afq.files.fibers.segmented{varargin{1}}) || ...
            ~ischar(afq.files.fibers.segmented{varargin{1}});
    case{'cleanfibers' 'cleanfibergroups' 'cleanfg'}
        val = logical(afq.overwrite.fibers.clean(varargin{1})) || ...
            isempty(afq.files.fibers.clean{varargin{1}}) || ...
            ~ischar(afq.files.fibers.clean{varargin{1}});
    case{'computevals' 'computeprofiles' 'computetractprofiles' 'compute'}
        % User wants to overwrite values for this subject
        val = logical(afq.overwrite.vals(varargin{1})) || ...
            % No values have been computed yet
            isempty(afq.vals.fa) || ...
            % Values have not yet been computed for this subject
            size(afq.vals.fa{1},1) < varargin{1};
            
    otherwise
        error('Uknown afq parameter');
end