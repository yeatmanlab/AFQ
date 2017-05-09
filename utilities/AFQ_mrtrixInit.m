function files = AFQ_mrtrixInit(dt6, ...
                                lmax, ...
                                mrtrix_folder, ...
                                mrtrixVersion, ...
                                multishell, ...
                                tool)
% function files = AFQ_mrtrixInit(dt6, lmax, mrtrix_folder)
% 
% Initialize an mrtrix CSD analysis
%
% This fucntion computes all the files needed to use mrtrix_track. 
%
% Parameters
% ----------
% dt6: string, full-path to an mrInit-generated dt6 file.
% T1nii: path to the acpc-ed T1w nii used at the beginning. 
% lmax: The maximal harmonic order to fit in the spherical deconvolution (d
%       model. Must be an even integer. This input determines the
%       flexibility  of the resulting model fit (higher values correspond
%       to more flexible models), but also determines the number of
%       parameters that need to be fit. The number of dw directions
%       acquired should be larger than the number of parameters required.
% mrtrix_folder; Name of the output folder
%
% Notes
% -----
% This performs the following operations:
%
% 1. Convert the raw dwi file into .mif format
% 2. Convert the bvecs, bvals into .b format
% 3. Convert the brain-mask to .mif format 
% 4. Fit DTI and calculate FA and EV images
% 5. Estimate the response function for single fibers, based on voxels with
%    FA > 0.7
% 6. Fit the CSD model. 
% 7. Convert the white-matter mask to .mif format. 
% 
% For details: 
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
% 


if notDefined('mrtrix_folder'), mrtrix_folder = 'mrtrix'; end
if notDefined('lmax'), lmax = 8; end
% Loading the dt file containing all the paths to the fiels we need.
dt_info = load(dt6);

% Check if this is correct, dt_info.files has some relative and absolute
% paths, it doesn't make sense. 
% I cannot recover the position of my original t1 from this information, so I
% copied it to the 'session' folder, in SubName/dmri 
%                 b0: 'dti90trilin/bin/b0.nii.gz'
%          brainMask: 'dti90trilin/bin/brainMask.nii.gz'
%             wmMask: 'dti90trilin/bin/wmMask.nii.gz'
%            tensors: 'dti90trilin/bin/tensors.nii.gz'
%                gof: 'dti90trilin/bin/gof.nii.gz'
%           outliers: 'dti90trilin/bin/outliers.nii.gz'
%                 t1: 't1_std_acpc.nii.gz'
%       alignedDwRaw: '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/data_aligned_...'
%     alignedDwBvecs: '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/data_aligned_...'
%     alignedDwBvals: '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/data_aligned_...'


% Note GLU: this code is assuming there is a 'raw' folder. In my case there
% is, but only with the original .nii-s converted from dicoms and the bvecs
% and bvals. The t1-s are in other path with the rest of the anat files.
% Furthermore, the assumption that the 'raw' file is above the dt6 filename
% breaks the code as it is duplicating the whole pathnames. 
% Example: mrtrix_dir = /bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002//bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/dti90trilin/mrtrixi/
% I fixed this (anf the initial part of the name issue as well)
% I still don't understand the use case. I understand that the mrtrix
% folder should be at the same level as the dt6.mat file, which defines
% every subject analysis, so using the bde_ code, it should be below the
% dti90trilin folder. I understand that the piece of filename that wants to
% be saved is the 'data_aligned_trilin_noMEC' part.


% Strip the file names out of the dt6 strings. 
% dwRawFile = dt_info.files.alignedDwRaw;
dwRawFile = fullfile(dt_info.params.rawDataDir, strcat(dt_info.params.rawDataFile,'.gz'));


% This line removes the extension of the file (.nii.gz) and mainaints de path
fname_trunk = dwRawFile(1:strfind(dwRawFile,'.')-1);
% With this code we can separate the rest
[pathDwRawFile, fnameDwRawFile] = fileparts(fname_trunk);

% In the mrtrix_folder argument we already have the path to the mrtrix
% folder
mrtrix_dir = mrtrix_folder;

% Assuming in 'session' we want the subject_name/dmri64 or whatever
% session = pathDwRawFile;
session = dt_info.params.rawDataDir;
% And in fname_trunk we want the whole path and the beginning of the
% filename
fname_trunk = [mrtrix_folder filesep fnameDwRawFile]; 

if ~exist(mrtrix_dir, 'dir')
    mkdir(mrtrix_dir)
end

% Build the mrtrix file names.
files = AFQ_mrtrix_build_files(fname_trunk, lmax, multishell);

% Check wich processes were already computed and which ons need to be doen.
computed = mrtrix_check_processes(files);

% Convert the raw dwi data to the mrtrix format: 
if (~computed.('dwi'))
    AFQ_mrtrix_mrconvert(dwRawFile, ...
                         files.dwi, ...
                         0, ...
                         0, ...
                         mrtrixVersion); 
end


% This file contains both bvecs and bvals, as per convention of mrtrix
% if (~computed.('b'))
    bvecs = fullfile(dt_info.params.rawDataDir, strcat(fnameDwRawFile, '.bvecs'));
    bvals = fullfile(dt_info.params.rawDataDir, strcat(fnameDwRawFile, '.bvals'));
%   bvecs = dt_info.files.alignedDwBvecs;
%   bvals = dt_info.files.alignedDwBvals;
    mrtrix_bfileFromBvecs(bvecs, bvals, files.b);
% end

% Convert the brain mask from mrDiffusion into a .mif file: 
if (~computed.('brainmask'))
  brainMaskFile = fullfile(session, dt_info.files.brainMask); 
  AFQ_mrtrix_mrconvert(brainMaskFile, ...
                       files.brainmask, ...
                       false, ...
                       false, ...
                       mrtrixVersion); 
end

% Generate diffusion tensors:
if (~computed.('dt'))
  AFQ_mrtrix_dwi2tensor(files.dwi, ...
                        files.dt, ...
                        files.b,...
                        0, ...
                        mrtrixVersion);
end

% Get the FA from the diffusion tensor estimates: 
if (~computed.('fa'))
  AFQ_mrtrix_tensor2FA(files.dt, ...
                       files.fa, ...
                       files.brainmask, ...
                       0, ...
                       mrtrixVersion);
end

% Generate the eigenvectors, weighted by FA: 
if  (~computed.('ev'))
  AFQ_mrtrix_tensor2vector(files.dt, files.ev, files.fa,0,mrtrixVersion);
end

% Estimate the response function of single fibers: 
if (~computed.('response'))
  AFQ_mrtrix_response(files.brainmask, ...
                      files.fa, ...
                      files.sf, ...
                      files.dwi,...
                      files.response, ...
                      files.b, ...
                      [], ... %this is the threshold string, it was missing!
                      false, ... % this is show figure
                      [], ... %bckground
                      8, ... %lmax
                      false, ... %verbose
                      mrtrixVersion) 
end

% Create a white-matter mask, tracktography will act only in here.
if (~computed.('wmMask'))
  wmMaskFile = fullfile(session, dt_info.files.wmMask);
  AFQ_mrtrix_mrconvert(wmMaskFile, ...
                       files.wmMask, ...
                       [], ...
                       0, ...
                       mrtrixVersion)
end



if ~multishell
    % Compute the CSD estimates: 
    if (~computed.('csd'))  
      disp('The following step takes a while (a few hours)');                                  
      AFQ_mrtrix_csdeconv(files.dwi, ...
                          files.response, ...
                          lmax, ...
                          files.csd, ... %out
                          files.b, ... %grad
                          files.brainmask,... %mask
                          false,... % Verbose
                          mrtrixVersion)
    end

else
    % Create the 5tt file from the same ac-pc-ed T1 nii we used in the other steps: 
    if (~computed.('tt5')) && (mrtrixVersion > 2)
        inputFile = [];
        if strcmp(tool, 'fsl')
            inputFile = fullfile(session, dt_info.files.t1);
            if ~(exist(inputFile, 'file') == 2)
                error(['Cannot find T1, please copy it to ' session]);
            end    
        % Find and aseg file and better if it is an aparc+aseg one. 
        % Select the first aparc if there are several *aseg* files.
        % It can take mgz or nii or mif
        else
           asegFiles = dir(fullfile(session,'*aseg*'));
           for ii = 1:length(asegFiles)
               if length(strfind(asegFiles(ii).name, 'aseg')) > 0
                   inputFile = fullfile(session, asegFiles(ii).name);
               end
               if length(strfind(asegFiles(ii).name, 'aparc')) > 0
                   inputFile = fullfile(session, asegFiles(ii).name);
               end
           end
            if ~(exist(inputFile, 'file') == 2)
                disp(['inputFile = ' inputFile]); 
                error(['Cannot find aseg file, please copy it to ' session]);
            end     
        end        

        % TODO: update directory structure to point to FS files.

        AFQ_mrtrix_5ttgen(inputFile, ...
                          files.tt5, ...
                          0, ...
                          0, ...
                          mrtrixVersion,...
                          tool);
    end
    
    % Create per tissue response function estimation
    % Not using the other response function, we already have the masks
    if (~computed.('wmResponse')) && (mrtrixVersion > 2)
        cmd_str = ['dwi2response msmt_5tt ' ...
                    files.dwi ' ' files.tt5 ' ' ...
                    files.wmResponse ' ' files.gmResponse ' ' files.csfResponse ' ' ...
                    '-grad ' files.b];
                    
        AFQ_mrtrix_cmd(cmd_str, 0, 0,mrtrixVersion);
    end
    
    
    % Compute the CSD estimates: 
    if (~computed.('csd'))   && (mrtrixVersion > 2)
      disp('The following step takes a while (a few hours)');                                  
      AFQ_mrtrix_csdeconv_msmt(files.dwi, ...
                              files.wmResponse, ...
                              files.gmResponse, ...
                              files.csfResponse, ...
                              lmax, ...
                              files.wmCsd, ...
                              files.gmCsd, ...
                              files.csfCsd, ...
                              files.b, ...
                              files.brainmask, ...
                              0, ...
                              0, ...
                              mrtrixVersion)
    end
    
    % RGB tissue signal contribution maps
    if (~computed.('vf'))  && (mrtrixVersion > 2)
         % mrconvert -coord 3 0 wm.mif - | mrcat csf.mif gm.mif - vf.mif
        cmd_str = ['mrconvert -coord 3 0 ' files.wmCsd ' - | ' ...
                   'mrcat ' files.gmCsd ' ' files.csfCsd ' - ' files.vf];           
        AFQ_mrtrix_cmd(cmd_str, 0, 0,mrtrixVersion);
    end
    
end



