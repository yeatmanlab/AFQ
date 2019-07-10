function files = AFQ_mrtrix_build_files(fname_trunk,lmax,compute5tt, multishell)
% Builds a structure with the names of the files that the MRtrix commands
% will generate and need.
%
% files = mrtrix_build_files(fname_trunk,lmax)
%
% Franco Pestilli, Ariel Rokem, Bob Dougherty Stanford University
% GLU July.2016 added the T1 and tt5 file types


% Convert the raw dwi data to the mrtrix format: 
files.dwi = strcat(fname_trunk,'_dwi.mif');

% This file contains both bvecs and bvals, as per convention of mrtrix
files.b     = strcat(fname_trunk, '.b');

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
files.wmMask    = strcat(fname_trunk, '_wmMask.mif');

% Compute the CSD estimates: 
files.csd = strcat(fname_trunk, sprintf('_csd_lmax%i.mif',lmax)); 

% Create tissue type segmentation to be used in multishell or ACT

if compute5tt>0 || mulitshell>0
files.tt5 = strcat(fname_trunk, '_5tt.mif');
files.gmwmi = strcat(fname_trunk, '_5tt_gmwmi.mif');
end

if multishell>0
    % Create per tissue type response file
    files.wmResponse = strcat(fname_trunk, '_wmResponse.txt');
    files.gmResponse = strcat(fname_trunk, '_gmResponse.txt');
    files.csfResponse = strcat(fname_trunk, '_csfResponse.txt');
    % Compute the CSD estimates: 
    files.wmCsd  = strcat(fname_trunk, sprintf('_wmCsd_lmax%i.mif',lmax)); 
    files.gmCsd  = strcat(fname_trunk, sprintf('_gmCsd_lmax%i.mif',lmax)); 
    files.csfCsd = strcat(fname_trunk, sprintf('_csfCsd_lmax%i.mif',lmax)); 
    % RGB tissue signal contribution maps
    files.vf = strcat(fname_trunk, '_vf.mif');
end
end