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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% For now don't run mrQ
if notDefined('run_mrQ')
    run_mrQ = 1;
end
% But do run AFQ
if notDefined('run_AFQ')
    run_AFQ = 1;
end 

% Check if dataDirs and outdirs are cell arrays and if not make it so
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
    try
        %% Run mrQ
        if run_mrQ == 1
            
            % Create the initial structure
            outDir = fullfile(outList{ii},'mrQ');
            
            if ~exist(outDir,'dir'); mkdir(outDir); end
            mrQ = mrQ_Create(inList{ii},[],outDir);
            
            % Set other parameters
            mrQ = mrQ_Set(mrQ,'sub',num2str(ii));
            mrQ = mrQ_Set(mrQ,'proclus',0);
            mrQ = mrQ_Set(mrQ,'sungrid',1);
            mrQ = mrQ_Set(mrQ,'fieldstrength',3);
            
            % New input to automatically acpc align
            mrQ = mrQ_Set(mrQ,'autoacpc',1);
            
            % Specific arrange function for nimsfs
            mrQ = mrQ_arrangeData_nimsfs(mrQ);
            
            % RUN IT
            mrQ_run(mrQ.name);
            
            % Reload the new mrQ structure
            mrQ = mrQ_Create(outDir);
            
            % Get a path to the synthesized t1file from mrQ
            t1File = fullfile(mrQ.spgr_initDir,'T1wfs_4.nii.gz');
            
        end
        
        if run_AFQ == 1
            %% Run dtiInit for each dwi data set
            
            % Get a list of directories for the diffusion data
            dwiFiles = getDwiFilesStruct(inList{ii});
            
            % VERY IMPORTANT - this simply chooses the last series number!!!
            % See getDwiFilesStruct for more details.
            dwiData  = dwiFiles{numel(dwiFiles)};
            
            % dwiData might actually be a structure with multiple diffusion
            % files. See above (important).
            for jj = 1:numel(dwiData)
                
                % If no t1 image was provided then use a template
                if notDefined('t1File')
                    t1File = fullfile(AFQ_directories,...
                        'templates',...
                        'mni_icbm152_nlin_asym_09a_nifti',...
                        'mni_icbm152_t2_tal_nlin_asym_09a.nii');
                end
                
                % These are the parameters we used for the weston havens data
                params = dtiInitParams;
                params.eddyCorrect = 0;
                params.noiseCalcMethod = 'b0';
                params.numBootStrapSamples = 0;
                params.clobber = -1;
                params.fitMethod = 'ls';
                
                % This is where we'll save all the outputs
                params.dt6BaseName = fullfile(outList{ii},['dt' num2str(dwiData.directions)]);
                
                % Get a path to the dwi file.
                dwiFile          = dwiData.nifti;
                params.bvecsFile = dwiData.bvec;
                params.bvalsFile = dwiData.bval;
                
                % Run dtiInit
                dtPath{jj} = dtiInit(dwiFile,t1File,params);
                
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
                'sub_group',sub_group,'clip2rois',0,'outdir',outList{ii});
            
            
            if run_mrQ == 1
                
                % First we need to align the maps to the diffusion data.
                % The first map in the cell array is the t1 - which is the
                % first input argument to the register funciton, so we skip
                % that one here.
                maps = fieldnames(mrQ.maps); maps = maps(2:end);
                otherMapsPath = cell(1,numel(maps));
                for i = 1:numel(maps)
                    otherMapsPath{i} = mrQ.maps.(maps{i});
                end
                
                t1 = mrQ.maps.T1path;
                b0 = mrvFindFile('b0.nii.gz',outList{ii});
                %mapsdir = fulfile(params.dt6BaseName,'bin');
                mapsdir = mrvDirup(b0);
                
                % Perform the alignment of the maps to the dti data (b0)
                alignedMaps = mrQ_registerMap2DTI(b0,t1,otherMapsPath,mapsdir);
                
                % Get and set the fieldnames for those maps we want to analyze
                % with AFQ.
                % *** These need to be the _DTI aligned maps.
                % mapsNames = fieldnames(mrQ.maps);
                for jj = 1:numel(alignedMaps)
                    % afq = AFQ_set(afq, 'images', mrQ.maps.(mapsNames{jj}));
                    % This needs to be a cell array, or AFQ_set complains
                    image{1} = alignedMaps{jj};
                    afq = AFQ_set(afq, 'images', image);
                end
                
            end
            
            % Run AFQ
            afq = AFQ_run(sub_dirs, sub_group, afq);
            
            % Setup the valnames.
            valnames = fieldnames(afq.vals);
            
            % %%%%% remove those values we're not interested in:
            remlist = {'cl','curvature','torsion','volume','WF_map_2DTI','cT1SIR_map_2DTI'};
            for r = 1:numel(remlist)
                ind = cellfind(valnames,remlist{r});
                if ~isempty(ind)
                    valnames(ind) = [];
                end
            end
            
            % Load up saved controls data
            % afq_controls = load('/biac4/wandell/data/WH/analysis/AFQ_WestonHavens_Full.mat');
            afq_controls = load('/biac4/wandell/biac2/wandell2/data/WH/analysis/WH_database_current.mat');
            
            % Save out images for the webpage
            
            % Create an output directory of figures
            figsDir = fullfile(outList{ii},'figures');
            if ~exist(figsDir,'dir')
                mkdir(figsDir);
            end
            
            
            AFQ_PlotPatientMeans(afq, afq_controls.afq, valnames, 21:80, figsDir);
            
        end
        
    catch theCatch
        
        % If anything in the above code crashes, then let's note that there was
        % an error and write that error out to the .error file.
        fprintf('ERROR in %s\n\n', inList{ii});
        if exist(fullfile(inList{ii},'.processing'),'file')
            movefile(fullfile(inList{ii},'.processing'),fullfile(inList{ii},'.error'));
            fprintf('\n===== Error Message ====\n');
            fprintf('%s\n\n',theCatch.message);
            
            % Write it out to the file
            fid = fopen(fullfile(inList{ii},'.error'),'a');
            fprintf(fid,'\n===== Error Message ====\n');
            fprintf(fid,'%s\n',theCatch.message);
            fclose(fid);
        else
            fprintf('No Error file found!\n===== Error Message ====\n');
            fprintf('%s\n\n',theCatch.message);
        end
    end
end
