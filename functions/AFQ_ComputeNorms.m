function [norms, patient_data, control_data, afq] = AFQ_ComputeNorms(varargin)
% Compute tract diffusion profile norms from control group data
%
% [norms, patient_data, control_data, afq] = AFQ_ComputeNorms(afq)
% or
% [norms, patient_data, control_data] = AFQ_ComputeNorms(groupFA, groupMD, groupRD, groupAD, groupCL, sub_group)
%
% This function can be called in 2 ways. Either with an afq structure as
% the input or with cell arrays of all the values and a subject grouping
% variable. Data will be separated into patient and control data, norms
% will be calculated and returned and if an afq structure was passed in
% then it will be modified to contain the new information.
%
% Inputs:
% afq       = An afq structure containing values and tract profiles for all
%             subjects. See AFQ_Create and AFQ_run
% groupFA   = Cell array of FA values for each subject's fiber tracts
% groupRD   = Cell array of RD values for each subject's fiber tracts
% groupAD   = Cell array of AD values for each subject's fiber tracts
% groupCL   = Cell array of CL values for each subject's fiber tracts
% sub_group = 0 or 1 denoting whether subject's are patients (1) or
%             controls (0)
%
% Outputs:
% norms        = means and standard deviations for each point on each fiber
%                tract
% patient_data = Data for the patients
% control_data = Data for the controls
% afq          = afq structure with new fields containing the norms,
%                patient and control group data

if nargin == 1
    afq = varargin{1};
    % Compute the norms and add them to the afq structure. This call also
    % separates the patient and control data and adds these fields to the
    % afq structure
    afq = AFQ_set(afq,'norms');
    % Get the norms, patient and control data from the afq structure
    norms = AFQ_get(afq,'norms');
    patient_data = AFQ_get(afq,'patient_data');
    control_data = AFQ_get(afq,'control_data');
elseif nargin == 6
    % If no afq structure was provided than assume inputs are Tract Diffusion
    % Profiles
    groupFA = varargin{1}; groupMD = varargin{2}; groupRD = varargin{3};
    groupAD = varargin{4}; groupCL = varargin{5}; sub_group = varargin{6};
    
    % Loop over fiber tracts and compute norms
    for jj=1:length(groupFA)
        
        norms.meanFA(:,jj)  = nanmean(groupFA{jj}(sub_group==0,:));
        norms.sdFA(:,jj)    = nanstd(groupFA{jj}(sub_group==0,:));
        patient_data(jj).FA = groupFA{jj}(sub_group==1,:);
        control_data(jj).FA = groupFA{jj}(sub_group==0,:);
        
        norms.meanMD(:,jj)  = nanmean(groupMD{jj}(sub_group==0,:));
        norms.sdMD(:,jj)    = nanstd(groupMD{jj}(sub_group==0,:));
        patient_data(jj).MD = groupMD{jj}(sub_group==1,:);
        control_data(jj).MD = groupMD{jj}(sub_group==0,:);
        
        norms.meanRD(:,jj)  = nanmean(groupRD{jj}(sub_group==0,:));
        norms.sdRD(:,jj)    = nanstd(groupRD{jj}(sub_group==0,:));
        patient_data(jj).RD = groupRD{jj}(sub_group==1,:);
        control_data(jj).RD = groupRD{jj}(sub_group==0,:);
        
        norms.meanAD(:,jj)  = nanmean(groupAD{jj}(sub_group==0,:));
        norms.sdAD(:,jj)    = nanstd(groupAD{jj}(sub_group==0,:));
        patient_data(jj).AD = groupAD{jj}(sub_group==1,:);
        control_data(jj).AD = groupAD{jj}(sub_group==0,:);
        
        norms.meanCL(:,jj)  = nanmean(groupCL{jj}(sub_group==0,:));
        norms.sdCL(:,jj)    = nanstd(groupCL{jj}(sub_group==0,:));
        patient_data(jj).CL = groupCL{jj}(sub_group==1,:);
        control_data(jj).CL = groupCL{jj}(sub_group==0,:);
        
    end
else
    error('impropper number of input arguments');
end