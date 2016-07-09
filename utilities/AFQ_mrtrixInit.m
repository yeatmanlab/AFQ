function files = AFQ_mrtrixInit(dt6, T1nii, lmax, mrtrix_folder, mrtrixVersion)
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
% Edit GLU July 2016. 
% 1: I think that due to the 2 dir structure for the multishell
% DWI, the directory structure is broken. I am going to give every function the
% path it needs, right now it is not working. I will just comment the old
% version if I need it afterwards. 
% 2: 
% TODO: Make it multishell within mrTrix. http://mrtrix.readthedocs.io/en/latest/workflows/multi_tissue_csd.html
% 1.- Obtain the 5 tissue type-s segmentation (5tt)
% 2.- 
%
%
% Example data: 
% dt6 = '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/dti90trilin/dt6.mat'
% lmax = 4;
% mrtrix_folder = '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/dti90trilin/mrtrix'
% mrtrixVersion = 3;




afq.software.mrtrixVersion = mrtrixVersion;

% if notDefined('mrtrix_folder'), 
%     mrtrix_folder = 'mrtrix'; 
% end
if notDefined('lmax'), lmax = 8; end
% Loading the dt file containing all the paths to the fiels we need.
dt_info = load(dt6);

% Check if this is correct, dt_info.files has some relative and absolute
% paths, it doesn't make sense. 
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
% I am trying to
% fix this, the fix I did in mac would leave the mrtrix names withouth the
% initial part of the filename. 
% I am going to comment the code in order to be easier for you to revise
% the code. I still don't know the use case. I understand that the mrtrix
% folder should be at the same level as the dt6.mat file, which defines
% every subject analysis, so using the bde_ code, it should be below the
% dti90trilin folder. I understand that the piece of filename that wants to
% be saved is the 'data_aligned_trilin_noMEC' part.


% Strip the file names out of the dt6 strings. 
dwRawFile = dt_info.files.alignedDwRaw;
T1niiFile = dt_info.files.t1;


% This line removes the extension of the file (.nii.gz) and mainaints de path
fname_trunk = dwRawFile(1:strfind(dwRawFile,'.')-1);
% With this code we can separate the rest
[pathDwRawFile, fnameDwRawFile] = fileparts(fname_trunk);

% In the mrtrix_folder argument we already have the path to the mrtrix
% folder
mrtrix_dir = mrtrix_folder;

% Assuming in 'session' we want the subject_name/dmri64 or whatever
session = pathDwRawFile; 
% And in fname_trunk we want the whole path and the beginning of the
% filename
fname_trunk = [mrtrix_folder filesep fnameDwRawFile]; 



if ~exist(mrtrix_dir, 'dir')
    mkdir(mrtrix_dir)
end

% Build the mrtrix file names.
files = AFQ_mrtrix_build_files(fname_trunk,lmax);

% Check wich processes were already computed and which ons need to be doen.
computed = mrtrix_check_processes(files);

% Convert the raw dwi data to the mrtrix format: 
if (~computed.('dwi'))
    AFQ_mrtrix_mrconvert(dwRawFile, ...
                         files.dwi, ...
                         0, ...
                         0, ...
                         afq.software.mrtrixVersion); 
end

% Convert the acpc nii T1w data to the mrtrix format: 
if (~computed.('T1'))
    AFQ_mrtrix_mrconvert(fullfile(session, T1niiFile), ...
                         files.T1, ...
                         0, ...
                         0, ...
                         afq.software.mrtrixVersion); 
end


% Create the 5tt file from the T1.mif: 
if (~computed.('tt5')) && (afq.software.mrtrixVersion > 2)
    % do this: 5ttgen fsl/freesurfer T1.mif 5tt.mif
       AFQ_mrtrix_5ttgen(fullfile(session, T1niiFile), ...
                         files.tt5, ...
                         0, ...
                         0, ...
                         afq.software.mrtrixVersion,...
                         'freesurfer'); 
end

% This file contains both bvecs and bvals, as per convention of mrtrix
if (~computed.('b'))
  bvecs = dt_info.files.alignedDwBvecs;
  bvals = dt_info.files.alignedDwBvals;
  mrtrix_bfileFromBvecs(bvecs, bvals, files.b);
end

% Convert the brain mask from mrDiffusion into a .mif file: 
if (~computed.('brainmask'))
  brainMaskFile = fullfile(session, dt_info.files.brainMask); 
  AFQ_mrtrix_mrconvert(brainMaskFile, ...
                       files.brainmask, ...
                       false, ...
                       false, ...
                       afq.software.mrtrixVersion); 
end

% Generate diffusion tensors:
if (~computed.('dt'))
  AFQ_mrtrix_dwi2tensor(files.dwi, ...
                        files.dt, ...
                        files.b,...
                        0, ...
                        afq.software.mrtrixVersion);
end

% Get the FA from the diffusion tensor estimates: 
if (~computed.('fa'))
  AFQ_mrtrix_tensor2FA(files.dt, files.fa, files.brainmask,0,afq.software.mrtrixVersion);
end

% Generate the eigenvectors, weighted by FA: 
if  (~computed.('ev'))
  AFQ_mrtrix_tensor2vector(files.dt, files.ev, files.fa,0,afq.software.mrtrixVersion);
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
                      afq.software.mrtrixVersion) 
end

% Create a white-matter mask, tracktography will act only in here.
if (~computed.('wm'))
  wmMaskFile = fullfile(session, dt_info.files.wmMask);
  AFQ_mrtrix_mrconvert(wmMaskFile, ...
                       files.wm, ...
                       [], ...
                       0, ...
                       afq.software.mrtrixVersion)
end

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
                      afq.software.mrtrixVersion)
end




