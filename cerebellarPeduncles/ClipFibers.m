function [] = ClipFibers(afq,fgNames,roi1cNames,roi2cNames,tdir_clipped,type)

%% This function saves the clipped version of the cerebellar peduncles for visual inspection. 
%  Note: Clipping can be done using 'and' or 'not' methods. Clipped fibers will be saved as
%  separate files and not override the full fibers. 

%  Written by Maya Yablonski - June 2016 
%  Edited for publication - January 2019
    
%% Input variables:
%  afq          - An afq structure output from AFQ_run (see also AFQ_Create)
%  fgNames      - A cell array of the names of the cerebellar peduncles.
%                 The fiber groups (fg) should be in each subject's fibers directory (with the same name). 
%                 This can be either a .mat file or a .pdb file.
%  roi1cNames   - A cell array of the names of the first waypoint ROIs used for the cererbellar segmentation.
%                 These can either be .mat files that are in each subject's ROI directory 
%                 or they can be a nifti images containing the ROI defined in MNI space. 
%                 If the ROIs are nifti images then they will be warped into each subject's native space 
%                 and saved as a .mat file. This allows new segmentations to be added to AFQ
%                 simply by defining the ROIs on the MNI template
%  roi2cNames   - A cell array of the names of the second waypoint ROIs used for the cererbellar segmentation.
%                 Either .mat files that are in each subject's ROI directory or nifti images in MNI space)
%  tdir_clipped - The path to the cerebellar template ROIs
%  type         - Options: 'and' =
%                          'not' =

%% Clip cerebellar peduncles to ROIs used for segmentation and save

for ii = 1:length(roi1cNames)
    if strfind(roi1cNames{ii},'/') % if ROIs already contain full path
        roi1c{ii} = roi1cNames{ii};
        roi2c{ii} = roi2cNames{ii};
    else
        roi1c{ii} = fullfile(tdir_clipped,roi1cNames{ii});
        roi2c{ii} = fullfile(tdir_clipped,roi2cNames{ii});
    end
end

runsubs = AFQ_get(afq,'run subjects');
for ii = 1:length(fgNames)
    for jj= runsubs
        sdir = afq.sub_dirs(jj);
        cur_fibFile = fullfile(sdir,'fibers',fgNames(ii));
        if strfind(cur_fibFile{1},'.mat')
            cur_fibFile = [cur_fibFile{1}];
        else
            cur_fibFile = [cur_fibFile{1} '.mat'];
        end
        
        if exist(cur_fibFile,'file')
            cur_fg = dtiLoadFiberGroup(cur_fibFile);
            
            roiPath = fullfile(sdir, 'ROIs');
            roi1File = [roiPath{1} '/' roi1cNames{ii}(1:end-7) '.mat'];
            roi2File = [roiPath{1} '/' roi2cNames{ii}(1:end-7) '.mat'];
            
            % =============================================================================================== %
            % If ROIs were not transformed to native space before:
            
            if ~exist(roi1File,'file') || ~exist(roi2File,'file')
                % Path to the directory of the template used for transformation
                tdir = fullfile(fileparts(which('mrDiffusion.m')), 'templates');
                template = fullfile(tdir,'MNI_JHU_T2.nii.gz');
                % Load subject's dt6-file
                dtpath = AFQ_get(afq,'dt6path',jj); dt = dtiLoadDt6(dtpath);
                % Specify subject's directory
                subdir = fileparts(dtpath);
               % Check if there is a precomputed spatial normalization. If not, compute it.
                sn = AFQ_get(afq,'spatial normalization',jj);
                invDef = AFQ_get(afq,'inverse deformation',jj);
                if isempty(sn) || isempty(invDef)
                    alignIm = mrAnatHistogramClip(double(dt.b0),0.3,0.99);
                    [sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(alignIm, dt.xformToAcpc, template);
                end
                % Load template ROIs in MNI space and transform them into subject's native space
                [~, ~, roi1]=dtiCreateRoiFromMniNifti(dt.dataFile, roi1c{ii}(1:end-7), invDef, 0);
                [~, ~, roi2]=dtiCreateRoiFromMniNifti(dt.dataFile, roi2c{ii}(1:end-7), invDef, 0);
                % Save ROIs as .mat files. Don't overwrite if ROIs already exist in subject's directory
                if ~exist(roi1File,'file')
                    dtiWriteRoi(roi1,roi1File);
                end
                if ~exist(roi2File,'file')
                    dtiWriteRoi(roi2,roi2File);
                end
            else
                [~,roi1N,suf1] = fileparts(roi1File);
                [~,roi2N,suf2] = fileparts(roi2File);
                fprintf('\n %s and %s exist for subject %d\n',[roi1N suf1],[roi2N suf2],jj)
            end
           % =============================================================================================== % 
            
            roi1=dtiReadRoi(roi1File);
            roi2=dtiReadRoi(roi2File);
            
            clippedFileName = [cur_fibFile(1:end-4) '-clip' type '.mat'];
            if ~exist(clippedFileName,'file')
                if numel(cur_fg.fibers) > 0
                    if strcmp(type,'and')
                        fgClipped = cp_dtiClipFiberGroupToROIs(cur_fg,roi1,roi2,[], {type});
                        % Note: First, we use "cp_dtiClipFiberToROIs" to be consistent with the way AFQ computes tract-profiles.
                        % Then we add the NOT ROIs using "dtiIntersectFibersWithROI"
                    elseif strcmp(type,'not')
                        fgClipped = cp_dtiIntersectFibersWithRoi([],{type},[],roi1,cur_fg,'first');
                        fgClipped = cp_dtiIntersectFibersWithRoi([],{type},[],roi2,fgClipped,'last');
                    else
                        fprintf(1,'Unknown filetype, using AND as default.\n');
                        fgClipped = cp_dtiClipFiberGroupToROIs(cur_fg,roi1,roi2,[], {'and'});
                    end
                    
                    % Change name of fiber group (fg) in fg.name -
                    % otherwise, it may cause inconsistencies when adding cerbellar peduncles to afq structure
                    fgClipped.name = [fgNames{ii}(1:end-4) '-clip' type '.mat'];
                    dtiWriteFiberGroup(fgClipped,clippedFileName);
                    fprintf(1,['Clipping fiber Group ' cur_fibFile ', subject ' num2str(jj) '\n']);
                else
                    fprintf(1,['Not enough fibers in ' cur_fibFile ', subject ' num2str(jj) '\n']);
                end
            else
                fprintf(1,['Fiber Group ' cur_fibFile 'already exists for subject ' num2str(jj) '\n']);
            end
        else
            fprintf(1,['Fiber Group ' cur_fibFile 'not found for subject ' num2str(jj) '\n']);
        end
    end
end
end
