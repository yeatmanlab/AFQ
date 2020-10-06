function  [SuperFiber, fgResampled, TractProfile,  T1, TV]=babyFatTractQmr(dwiDir, sessid, runName, fgName)

% dwiTractQmri(dwiDir, sessid, runName, fgName)
% fgName: Name of fg file

fprintf('Compute tract Qmr for (%s,%s,%s)\n',sessid,runName,fgName);
runDir = fullfile(dwiDir,sessid,runName,'dti94trilin');
afqDir = fullfile(runDir, 'fibers','afq');
numNodes = 100; clip2rois = 0;
weighting='mean';
% Load fg array
fgFile = fullfile(afqDir, fgName);
if exist(fgFile,'file')
    load(fgFile,'fg');
    roifg=fg;
    
    vent=niftiRead('ventricles_edited.nii.gz')
    SEIR=niftiRead('SEIR_T1_acpc.nii.gz')
    
    SEIR_masked=SEIR
    SEIR_masked.data(vent.data>0)=nan;
    
    niftiWrite(SEIR_masked,'SEIR_T1_acpc_masked_ventr.nii.gz')
    clear('SEIR','vent')
   
    t1=SEIR_masked;
    
    %Check t1 header
    if ~all(t1.qto_xyz(:) == t1.sto_xyz(:))
        t1 = niftiCheckQto(t1);
    end
    
    roiDir=fullfile(runDir,'babyAFQROIs');
    % Compute a Tract t1
    [T1, md, rd, ad, cl, volume, TractProfile, SuperFiber, fgResampled]  = babyAFQ_ComputeTractPropertiesNanMean(roifg,t1, numNodes, clip2rois, runDir,roiDir,[],weighting,[]);
                                                                           %babyAFQ_ComputeTractProperties(fg_classified,dt,numberOfNodes,clip2rois,subDir, roiDir, fgnum, weighting, afq)
    %AFQ_ComputeTractProperties(fg_classified,dt,[numberOfNodes=30],[clip2rois=1],[subDir], [weighting = 1], afq)
    % Compute diffusion properties along the trajectory of the fiber groups.
    % Read the tv file
    vent=niftiRead('ventricles_edited.nii.gz')
    cmd_str=['mrconvert -force ' fullfile(runDir,'mrtrix/dwi_processed_aligned_trilin_noMEC_md_dki.mif') ' ' fullfile(runDir,'mrtrix/md_dki.nii.gz')]
    system(cmd_str);
    
    mdFile = fullfile(runDir,'mrtrix/md_dki.nii.gz');
    cmd_str=['mri_convert ' mdFile ' ' mdFile ' --reslice_like ventricles_edited.nii.gz']
    system(cmd_str)
    md = niftiRead(mdFile);

    md_masked=md;
    md_masked.data(vent.data>0)=nan;
    
    niftiWrite(md_masked,'md_acpc_masked_ventr.nii.gz')
    clear('md','vent')
    
    md=md_masked;
    % Check tv header
    if ~all(md.qto_xyz(:) == md.sto_xyz(:))
        md = niftiCheckQto(md);
    end
    
    % Compute a Tract tv
     [MD, md, rd, ad, cl, volume, TractProfile, SuperFiber, fgResampled] = babyAFQ_ComputeTractPropertiesNanMean(roifg,md, numNodes, clip2rois, runDir,roiDir,[],weighting,[]);
       
     
    cmd_str=['mrconvert -force ' fullfile(runDir,'mrtrix/dwi_processed_aligned_trilin_noMEC_fa_dki.mif') ' ' fullfile(runDir,'mrtrix/fa_dki.nii.gz')]
    system(cmd_str); 
    vent=niftiRead('ventricles_edited.nii.gz')
        
    faFile = fullfile(runDir,'mrtrix/fa_dki.nii.gz');
    cmd_str=['mri_convert ' faFile ' ' faFile ' --reslice_like ventricles_edited.nii.gz']
    system(cmd_str)
    
    fa = niftiRead(faFile);
    
    fa_masked=fa;
    fa_masked.data(vent.data>0)=nan;
    
    niftiWrite(fa_masked,'fa_acpc_masked_ventr.nii.gz')
    clear('fa','vent')
    
    fa=fa_masked;
    % Check tv header
    if ~all(fa.qto_xyz(:) == fa.sto_xyz(:))
        fa = niftiCheckQto(fa);
    end
    % Compute a Tract tv
     [FA, md, rd, ad, cl, volume, TractProfile, SuperFiber, fgResampled] = babyAFQ_ComputeTractPropertiesNanMean(roifg,fa, numNodes, clip2rois, runDir,roiDir,[],weighting,[]);
       
     
%      if size(T1,2)>1
%        T1=T1(:,num)
%        SuperFiber=SuperFiber(num);
%        TractProfile=TractProfile(num);
%    end
    % save tract Qmr
    T1AcrNodes=T1;
    MdAcrNodes=MD;
    FaAcrNodes=FA;
    
    T1Avg=nanmean(T1);
    MdAvg=nanmean(MD);
    FaAvg=nanmean(FA);
    
    tractFile = fullfile(afqDir, ['TractQmr_masked_ventr_', fgName]);
    save(tractFile,'T1AcrNodes','T1Avg','MdAcrNodes','MdAvg','FaAcrNodes','FaAvg');
end

            