function [norms patient_data control_data] = AFQ_ComputeNorms(groupFA, groupMD, groupRD, groupAD, groupCL, sub_group)
%Compute tract diffusion profile norms from control group data
%
% [norms patient_data control_data] = AFQ_ComputeNorms(groupFA, groupMD, groupRD, groupAD, groupCL, sub_group)
%
% Inputs:
% groupFA   = Cell array of FA values for each subject's fiber tracts
% groupRD   = Cell array of RD values for each subject's fiber tracts
% groupAD   = Cell array of AD values for each subject's fiber tracts
% groupCL   = Cell array of CL values for each subject's fiber tracts
% sub_group = 0 or 1 denoting whether subject's are patients (1) or
%             controls (0)
%
% Outputs:
% norms     = means and standard deviations for each point on each fiber
%             tract
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