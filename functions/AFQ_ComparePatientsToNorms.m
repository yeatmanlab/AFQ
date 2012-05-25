function [abn abnTracts] = AFQ_ComparePatientsToNorms(patient_data, norms, cutoff, property, comp)
% Determine which premies fall outside the defined percentile band for each
% tract
%
% [abn abnTracts] = AFQ_ComparePatientsToNorms(patient_data, norms, cutoff, property, comp)
%
% Inputs:
% patient_data = Fiber tract diffusion profiles for the patients.  See
%                AFQ_ComputeTractProperties and AFQ_ComputeNorms
% norms        = Normative fiber tract diffusion profiles see AFQ_ComputeNorms
% cutoff       = Percentile bands that define the range to be considered
%                normal
% property     = String defining what property to analyze.  For example
%                'FA' or 'RD'
% comp         = Defines whether to compare tract means (comp = 'mean') or
%                compare pointwise along the tract profile
%                (comp ='profile').  The default is to compare along the
%                tract profiles.
%
%
% Outputs:
% abn          = A 1 x N vector where N is the number of patients.
%                Each patient that is abnormal on at least one tract is
%                marked with a 1 and each subject that is normal on every
%                tract is marked with a 0. The criteria for abnormal is
%                defined in afq.params.cutoff.  See AFQ create
%
% abnTracts    = An M by N matrix where M is the number of subjects and N
%                is the number of tracts. Each row is a subject and each
%                column is a tract.  1 means that tract was abnormal for
%                that subject and 0 means it was normal.
%
% Copyright Jason D. Yeatman. Vista Team 2011
%

if ~exist('comp','var') || isempty(comp)
    comp = 'profile';
end
cutZ=norminv(cutoff*.01);
for ii=1:size(patient_data(1).FA,1)
    for jj=1:length(patient_data)
        
        % Collect the property of interest and the relevant norms
        switch(property)
            case {'FA' 'fa' 'fractional anisotropy'}
                vals     = patient_data(jj).FA(ii,:)';
                val_mean = norms.meanFA(:,jj);
                val_sd   = norms.sdFA(:,jj);
            case {'MD' 'md' 'mean diffusivity'}
                vals     = patient_data(jj).MD(ii,:)';
                val_mean = norms.meanMD(:,jj);
                val_sd   = norms.sdMD(:,jj);
            case {'RD' 'rd' 'radial diffusivity'}
                vals     = patient_data(jj).RD(ii,:)';
                val_mean = norms.meanRD(:,jj);
                val_sd   = norms.sdRD(:,jj);
            case {'AD' 'ad' 'axial diffusivity'}
                vals     = patient_data(jj).AD(ii,:)';
                val_mean = norms.meanAD(:,jj);
                val_sd   = norms.sdAD(:,jj);
        end
        
        % Each subject gets a 0 or 1 to signify if any value for each
        % tract falls outside the desired confidence interval
        if strcmpi(comp, 'mean')
            % Abnormality is defined based on tract mean
            abnTracts(ii,jj) = nanmean(vals) > mean(val_mean + max(cutZ) .* val_sd) || ...
                nanmean(vals) < mean(val_mean + min(cutZ) .* val_sd);
        elseif strcmpi(comp, 'profile')
            % Abnormality is defined as any point on the tract profile
            % passing outside the defined confidence interval
            locs = vals > val_mean + max(cutZ) .* val_sd | ...
                vals < val_mean + min(cutZ) .* val_sd;
            abnTracts(ii,jj) = sum(locs) > 0;
        end
    end
end
% Find subjects that are abnomral on at least 1 tract
abn=sum(abnTracts,2)>0;