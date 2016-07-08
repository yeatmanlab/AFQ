function files = AFQ_mrtrixInit(dt6, lmax, mrtrix_folder, mrtrixVersion)
% function files = AFQ_mrtrixInit(dt6, lmax, mrtrix_folder)
% 
% Initialize an mrtrix CSD analysis
%
% This fucntion computes all the files needed to use mrtrix_track. 
%
% Parameters
% ----------
% dt6: string, full-path to an mrInit-generated dt6 file. 
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



afq.software.mrtrixVersion = mrtrixVersion;

% if notDefined('mrtrix_folder'), 
%     mrtrix_folder = 'mrtrix'; 
% end
if notDefined('lmax'), lmax = 8; end
% Loading the dt file containing all the paths to the fiels we need.
dt_info = load(dt6);

% Strip the file names out of the dt6 strings. 
dwRawFile = dt_info.files.alignedDwRaw;


% fname_trunk = dwRawFile(1:strfind(dwRawFile,'.')-1);
% raw_idx = strfind(fname_trunk,'raw');
% if isempty(raw_idx)
%     % if the raw directory was not actually named raw then assume it is
%     % right above the file name
%     [~,rname]=fileparts(fileparts(fname_trunk));
%     raw_idx = strfind(fname_trunk,rname);
%     %error('Could not find raw directory')
% end
% session = fname_trunk(1:raw_idx-1);
% fname_trunk = [session, mrtrix_folder, fname_trunk(raw_idx+3:end)]; 
% file_sep_idx = strfind(fname_trunk, filesep);
% mrtrix_dir = fname_trunk(1:file_sep_idx(end)); 

% The previous code is not working for my case (GLU)
[fname_trunk, fname, fext] = fileparts(dwRawFile);
mrtrix_dir = mrtrix_folder;
session = fname_trunk;


if ~exist(mrtrix_dir, 'dir')
    mkdir(mrtrix_dir)
end

% Build the mrtrix file names.
files = mrtrix_build_files([mrtrix_dir filesep],lmax);

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

% Convert the raw dwi data to the mrtrix format: 
if (~computed.('T1'))
    AFQ_mrtrix_mrconvert(dwRawFile, ...
                         files.dwi, ...
                         0, ...
                         0, ...
                         afq.software.mrtrixVersion); 
end


% Convert the raw dwi data to the mrtrix format: 
if (~computed.('tt5'))
    AFQ_mrtrix_mrconvert(dwRawFile, ...
                         files.tt5, ...
                         0, ...
                         0, ...
                         afq.software.mrtrixVersion); 
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



%%%%%%%%%%%%%%%%%%%%%%%%%%
% mrtrix_check_processes %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function computed = mrtrix_check_processes(files)
%
% Check which mrtrix proceses were computed. Returns an array of 0's and
% 1's indicating which process was computed (1) and which one needs to be
% computed (0)
%

fields = fieldnames(files);
for ii = 1:length(fields)
  if exist(files.(fields{ii}),'file') == 2
    computed.(fields{ii}) = 1;
  else
    computed.(fields{ii}) = 0;
  end
end


%%%%%%%%%%%%%%%%%%%%%%
% mrtrix_build_files %
%%%%%%%%%%%%%%%%%%%%%%
function files = mrtrix_build_files(fname_trunk,lmax)
%
% Builds a structure with the mrtrix file names.
%

% Convert the raw dwi data to the mrtrix format: 
files.dwi = strcat(fname_trunk,'_dwi.mif');

% This file contains both bvecs and bvals, as per convention of mrtrix
files.b     = strcat(fname_trunk, '_b');

% Convert the brain mask from mrDiffusion into a .mif file: 
files.brainmask = strcat(fname_trunk,'_brainmask.mif');

% Generate diffusion tensors:
files.dt = strcat(fname_trunk, '_dt.mif');

% Get the FA from the diffusion tensor estimates: 
files.fa = strcat(fname_trunk, '_fa.mif');

% Generate the eigenvectors, weighted by FA: 
files.ev = strcat(fname_trunk, '_ev.mif');

% Estimate the response function of single fibers: 
files.sf = strcat(fname_trunk, '_sf.mif');
files.response = strcat(fname_trunk, '_response.txt');

% Create a white-matter mask, tracktography will act only in here.
files.wm    = strcat(fname_trunk, '_wm.mif');

% Compute the CSD estimates: 
files.csd = strcat(fname_trunk, sprintf('_csd_lmax%i.mif',lmax)); 

% Create a coregistered T1 data set from the same subject
files.T1    = strcat(fname_trunk, '_T1.mif');

% Create tissue type segmentation to be used in multishell: 
files.tt5 = strcat(fname_trunk, sprintf('_5tt.mif',lmax)); 





