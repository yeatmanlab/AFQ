function [status, results] = AFQ_mrtrix_tensor2FA(in_file, out_file, ...
                                         mask_file, verbose, mrtrixVersion)

%
% Calculate diffusion tensors.
%
% Parameters
% ----------
% in_file: The name of a dti file in .mif format
% out_file: The name of the resulting fa file in .mif format
% mask_file: The name of a mask file in .mif format (default to the entire
%             volume). 
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


if notDefined('verbose')
    verbose = true;
end

if notDefined('bkgrnd')
    bkgrnd = false;
end





if mrtrixVersion == 2
    func1Name = 'tensor2FA';
    func1NameOpt = '';
    func2Name = 'mrmult';
    func2NameOpt = '';
end
if mrtrixVersion == 3
    func1Name = 'tensor2metric';
    func1NameOpt = '-fa';
    func2Name = 'mrcalc';
    func2NameOpt = '-mult';
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