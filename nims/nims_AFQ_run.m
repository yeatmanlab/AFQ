function nims_AFQ_run(dataDirs, outdirs, run_mrQ, run_AFQ)
% Run dtiInit, AFQ and mrQ on raw data and save output images
%
% This function is intended to be launched when data is uploaded into the
% nims database (http://cni.stanford.edu/nims). Given a directory of nifti
% images created by nims, this function will run mrQ if the apropriate data
% is available, then it will run dtiInit.m if diffusion data is present and
% finally it will run AFQ on all the data that was input. Figures will be
% automatically saved out to be built into a webpage.
%
% Inputs:
%
% dataDirs - a cell array of subjects data directories to be run
% outdirs  - a cell array of directories to save each subject's processed
%            data and figures
% run_mrQ  - Logical denoting whether to run mrQ
% run_AFQ  - Logical denoting whether to run AFQ
%
% example:
%
% dataDirs = {'/scarlet/westonhavens/upload/testlab/20121011_1212_3425'};
% outdirs = {'/scarlet/westonhavens/results/testlab/20121011_1212_3425'};
% nims_AFQ_run(dataDirs, outdirs)

% For now don't run mrQ
if notDefined('run_mrQ')
    run_mrQ = 0;
end
% But do run AFQ
if notDefined('run_AFQ')
    run_AFQ = 1;
end
% Check if dataDirs and outdirs are cell arrays and if not make them cell
% arrays
if iscell(dataDirs)
    inList  = dataDirs;
else
    inList{1}  = dataDirs;
end
if iscell(outdirs)
    outList  = outdirs;
else
    outList{1}  = outdirs;
end

% Loop over all the input directories
for ii = 1:length(inList)
    
    %% This is all garbage that should be rewritten but will work as an example
    % List the files in the directory
    files = dir(inList{ii});
    % Check which file types are in each directory
    dwiData = {};
    spgrData = {};
    epiData = {};
    
    % Make a list of directories with dti data
    for jj = 1:length(files)
        if ~isempty(strfind(files(jj).name,'DTI_2mm_96dir'))
            dwiData = vertcat(dwiData,{fullfile(inList{ii},files(jj).name)});
        end
    end
    
    % Make a list of directories with SEIR data
    for jj = 1:length(files)
        if ~isempty(strfind(files(jj).name,'IR'))
            epiData = vertcat(epiData,{fullfile(inList{ii},files(jj).name)});
        end
    end
    
    % Make a list of directories with spgr data
    for jj = 1:length(files)
        if ~isempty(strfind(files(jj).name,'SPGR'))
            spgrData = vertcat(spgrData,{fullfile(inList{ii},files(jj).name)});
        end
    end
    
    %% Run mrQ
    if run_mrQ == 1
        mrQ=mrQ_Create(inList{ii});
        mrQ=mrQ_Set(mrQ,'SEIR',epiData);
        mrQ=mrQ_Set(mrQ,'SPGR',spgrData);
        mrQ=mrQ_Set(mrQ,'sub',num2str(ii));
        mrQ=mrQ_Set(mrQ,'proclus',0);
        mrQ=mrQ_Set(mrQ,'sungrid',0);
        % New input to automatically acpc align
        mrQ=mrQ_Set(mrQ,'autoacpc',1);
        mrQ_run(mrQ.name);
        % Reload the new mrQ structure
        mrQ=mrQ_Create(inList{ii});
        % Get a path to the synthasized t1file from mrQ
        t1File = fullfile(fileparts(mrQ.mapsDir,'T1wfs_4.nii.gz'));
    end
    
    if run_AFQ == 1
        %% Run dtiInit for each dwi data set
        for jj = 1:numel(dwiData)
            
            % If no t1 image was provided then use a template
            if notDefined('t1File')
                t1File = fullfile(AFQ_directories,'templates',...
                    'mni_icbm152_nlin_asym_09a_nifti','mni_icbm152_t2_tal_nlin_asym_09a.nii');
            end   
            
            % These are the parameters we used for the weston havens data
            params = dtiInitParams;
            params.eddyCorrect = 0;  params.noiseCalcMethod = 'b0';
            params.numBootStrapSamples = 0;
            params.clobber = -1;
            params.fitMethod = 'ls';
            % This is where we'll save all the outputs
            params.dt6BaseName = fullfile(outdirs{ii},'dt6');
            
            % Get a path to the dwi file. This is also a temporary call
            % that should be rewritten
            dwiFile = dir(fullfile(dwiData{jj},'*_1.nii.gz'));
            dwiFile = fullfile(dwiData{jj},dwiFile.name);
            % Get a path to the associated bvals and bvecs files
            params.bvecsFile = [prefix(prefix(dwiFile)) '.bvec'];
            params.bvalsFile = [prefix(prefix(dwiFile)) '.bval'];

            % Run dti init
            dtPath{jj}=dtiInit(dwiFile,t1File,params);
            
        end
        
        %% Run AFQ
        
        % Run afq for each dwi file that was input
        for jj = 1:numel(dtPath)
            sub_dirs{jj} = fileparts(dtPath{jj}{1});
        end
        % Just one group
        sub_group = ones(numel(sub_dirs),1);
        % Create afq structure
        afq = AFQ_Create('sub_dirs',sub_dirs,...
            'sub_group',sub_group,'clip2rois',0,'outdir',outdirs{ii});
        
        % Add in qMRI maps
        if run_mrQ == 1
            mapsNames = fieldnames(mrQ.maps);
            for jj = 1:numel(mapsNames)
                afq = AFQ_set(afq, 'images', mrQ.(mapsNames{jj}));
            end
        end
        
        % Run AFQ
        afq = AFQ_run(sub_dirs, sub_group, afq)
        
        %% Save out images for the webpage
        
        % Create an output directory of figures
        figsDir = fullfile(outdirs{ii},'figures');
        if ~exist(figsDir,'dir')
            mkdir(figsDir);
        end
        
        % Load up saved controls data
        afq_controls = load('/biac4/wandell/data/WH/analysis/AFQ_WestonHavens_Full.mat');
        
        % This is the function that makes figures.
        AFQ_PlotPatientMeans(afq, afq_controls.afq, [], 21:80, figsDir);
        
    end
end
