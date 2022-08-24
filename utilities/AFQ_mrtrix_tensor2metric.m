function [status, results] = AFQ_mrtrix_tensor2metric(in_file, out_file, ...
                                         mask_file, metric, verbose,mrtrixVersion)

%
% Calculate diffusion tensors.
%
% Parameters
% ----------
% in_file: The name of a dti file in .mif format
% out_file: The name of the resulting fa file in .mif format
% mask_file: The name of a mask file in .mif format (default to the entire
%             volume). 
% metric: metric to compute options are 'md', 'fa', 'ad', 'rd'
% note: only will work with mrTrix3
% 
% Returns
% -------
% status: whether (0) or not (non-zero) the operation succeeded
% results: the results of the operation in stdout
%
% Notes
% -----
% http://www.brain.org.au/software/mrtrix/tractography/preprocess.html
%
% Edited GLU 06.2016:
%        1.- Remove absolute paths
%        2.- Include mrTrix version
% 
% Edited ECK & XY 08.2022: 
%       1. only include functionality for mrTrix3
%       2. make more flexible for multiple metrics

if mrtrixVersion ~=3 
    error('mrTrix version must be 3.')
end 

if notDefined('verbose')
    verbose = true;
end

if notDefined('bkgrnd')
    bkgrnd = false;
end

func1Name = 'tensor2metric';
func1NameOpt = '-fa';
func2Name = 'mrcalc';
func2NameOpt = '-mult';

switch metric
    case 'fa'
        func1NameOpt = '-fa';
    case 'md'
        func1NameOpt = '-adc';
    case 'ad'
        func1NameOpt = '-ad';
    case 'rd'
        func1NameOpt = '-rd';
end

% If no mask file was provided, calculate this over the entire volume
if notDefined('mask_file')
    cmd_str = sprintf('%s %s %s %s', func1Name,    in_file, ... 
                                  func1NameOpt, out_file);
% Otherwise, use the mask file provided: 
else
    cmd_str = sprintf('%s %s %s - | %s - %s %s %s', ...
                       func1Name, in_file, func1NameOpt, ...
                       func2Name, mask_file, func2NameOpt, out_file);
end

% Send it to mrtrix:
[status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);