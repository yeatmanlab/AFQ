function afq = AFQ_ComputeSpatialNormalization(afq)
% Compute spatial normalization to MNI template for each subject
%
% afq = AFQ_ComputeSpatialNormalization(afq)
%

% Check which subjects should be run
runsubs = AFQ_get(afq,'run subjects');

if AFQ_get(afq, 'use ANTS')
    fprintf('\n Computing spatial normalization with ANTS\n')
    sub_dirs = AFQ_get(afq,'sub_dirs')
    for ii = runsubs
        % These will be the names of the output files
        warp = fullfile(sub_dirs{ii},'bin','b0Warp.nii.gz');
        invwarp = fullfile(sub_dirs{ii},'bin','b0InverseWarp.nii.gz');
        
        % Check if the already exist before runing ANTs
        if ~exist(warp,'file') || ~exist(invwarp,'file')
            b0File = fullfile(sub_dirs{ii},'bin','b0.nii.gz');
            tdir = fullfile(AFQ_directories,'templates','mni_icbm152_nlin_asym_09a_nifti');
            template = fullfile(tdir,'mni_icbm152_t2_tal_nlin_asym_09a.nii');
            ANTS_normalize(b0File,template);
        end
        
        % Set paths to these files in the afq structure
        afq = AFQ_set(afq,'ants warp', warp, ii);
        afq = AFQ_set(afq,'ants inverse warp', invwarp, ii);
    end
end