function afq = AFQ_AddFrontalAslantTract()
   
% This function allows the user to add bilateral Frontal Aslant Tracts (FAT) to an existing AFQ structure.
% Prior to the execution, the user should make sure that he knows where the AFQ structure is stored, and where the
% FAT ROIs are stored. 
% Note that the ROIs were defined on an MNI template template (ICBM 2009a Nonlinear Asymmetric template). 
% For details, see Kronfeld-Duenias et al., BSAF 2016; The full reference appears below.
%
% The Function was written by VKD on 20160303. 
% When used, please reference:
%'Kronfeld-Duenias, et al. "The frontal aslant tract underlies speech fluency in persistent developmental stuttering." 
% Brain Structure and Function 221.1 (2016): 365-381.‏'

    %% Load an AFQ structure and initiate the basice paprams
    % Choose an afq structure
    default_afq_file = '';
    [afqFileName,afqPathName] = uigetfile('*.mat','Select an AFQ file',default_afq_file);
    load(fullfile(afqPathName,afqFileName));

    if(~exist('afq'))
        msgbox('AFQ strucuture could not be succesfully loaded. Bailing out','Error');
        return;
    end
    
    % Should files be overriden?
    options.Interpreter = 'tex';
    options.Default = 'No';
    qstring = 'Overwrite existing files?';
    choice = questdlg(qstring,'Definitions','Yes','No',options)
    if(strcmp(choice,'Yes'))
        overwriteFiles = true;
    else
        overwriteFiles = false;
    end

    %% Add Left FAT
    % Choose a directory where the left FAT teplate ROIs are saved.
    tdir = uigetdir('','Choose the directory where the FAT ROIs are stored');

    roi1Name = fullfile(tdir,'l_aslant_roi1_x45.nii.gz');
    roi2Name = fullfile(tdir,'l_aslant_roi2_z45.nii.gz');
    fgName = 'L_ASLANT.mat';
    if ~exist(roi1Name,'file') || ~ exist(roi2Name,'file')
        msgbox('Left FAT ROIs could not be succesfully loaded. Bailing out','Error');
        return;
    end
    % Add the left FAT
    afq = AFQ_AddNewFiberGroup(afq,fgName,roi1Name,roi2Name,1,1,0,[],overwriteFiles);

    %% Add the right FAT
    % Choose a directory where the right FAT teplate ROIs are saved.
    roi1Name = fullfile(tdir,'r_aslant_roi1_x45.nii.gz');
    roi2Name = fullfile(tdir,'r_aslant_roi2_z45.nii.gz');
    fgName = 'R_ASLANT.mat';
    if ~exist(roi1Name,'file') || ~ exist(roi2Name,'file')
        msgbox('Right FAT ROIs could not be succesfully loaded. Bailing out','Error');
        return;
    end
    % Add the right FAT
    afq = AFQ_AddNewFiberGroup(afq,fgName,roi1Name,roi2Name,1,1,0,[],overwriteFiles);
    
    %% Save the updated AFQ
    % If requested, save the updated AFQ strucuture.
    options.Interpreter = 'tex';
    options.Default = 'Yes';
    qstring = 'Save the updated AFQ structure?';
    choice = questdlg(qstring,'Save','Yes','No',options)
    if(strcmp(choice,'Yes'))
        [file,path] = uiputfile('*.mat','Save AFQ structure');  
    end
    

end