function val = AFQ_get(afq, param, varargin)
% Get value from afq structure
%
% val = AFQ_get(afq, param, varargin)
%
% param list:
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
    otherwise
        error('Uknown afq parameter');
end