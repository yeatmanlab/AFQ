function [status, results] = AFQ_mrtrix_csdeconv(dwi_file, ...
                                                 response_file, ...
                                                 lmax, ...
                                                 out_file, ...
                                                 grad_file, ...
                                                 mask_file, ...
                                                 verbose,...
                                                 mrtrixVersion)
%function [status, results] = mrtrix_csdeconv(dwi_file, response_file, lmax, ...
%           out_file, grad_file, [mask_file=entire volume], [verbose=true])
%
% Fit the constrained spherical deconvolution model to dwi data 
%
% Parameters
% ----------
% dwi_file: The name of a dwi file in .mif format. 
% response_file: The name of a .txt fiber response function file. 
% lmax: The maximal harmonic order. 
% out_file: The resulting .mif file. 
% grad_file: a file containing gradient information in the mrtrix format. 
% mask_file: a .mif file containing a mask. Default: the entire volume. 
% verbose: Whether to display stdout to the command window. 
% 
% Returns
% -------
% status: whether (0) or not (1) the operation succeeded
% results: the results of the operation in stdout
%
% Notes
% -----
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
% 

if notDefined('verbose')
    verbose = true; 
end

if notDefined('bkgrnd')
    bkgrnd = false;
end


if mrtrixVersion == 2
    funcName = 'csdeconv';
end
if mrtrixVersion == 3
    funcName = 'dwi2fod csd';
end




if notDefined('mask_file')
    cmd_str = sprintf('%s %s %s -lmax %d %s -grad %s',...
                       funcName, dwi_file, response_file, ...
                       lmax, out_file, grad_file);

else 
    cmd_str = sprintf('%s %s %s -lmax %d -mask %s %s -grad %s',...
                       funcName, dwi_file, response_file, ...
                       lmax, ...
                       mask_file, out_file, ...
                       grad_file);
end 

[status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);



% Esto es lo que he lanzado en command line
% dwi2fod msmt_csd -mask data_aligned_trilin_noMEC_brainmask.mif data_aligned_trilin_noMEC_dwi.mif data_aligned_trilin_noMEC_wm.txt data_aligned_trilin_noMEC_wm.mif data_aligned_trilin_noMEC_gm.txt data_aligned_trilin_noMEC_gm.mif data_aligned_trilin_noMEC_csf.txt data_aligned_trilin_noMEC_csf.mif -grad data_aligned_trilin_noMEC.b




