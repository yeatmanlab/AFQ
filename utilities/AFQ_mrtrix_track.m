function [status, results, fg, pathstr] = AFQ_mrtrix_track(files, ...
                                                           roi, ...
                                                           mask, ...
                                                           mode, ...
                                                           nSeeds, ...
                                                           bkgrnd, ...
                                                           verbose, ...
                                                           clobber, ...
                                                           mrtrixVersion)
%
% function [status, results, fg, pathstr] = mrtrix_track(files, roi, mask, mode, nSeeds, bkgrnd, verbose)
%
% Provided a csd estimate, generate estimates of the fibers starting in roi 
% and terminating when they reach the boundary of mask
%
% Parameters
% ----------
% files: structure, containing the filenames generated using mrtrix_init.m
% roi: string, filename for a .mif format file containing the ROI in which
%      to place the seeds. Use the *_wm.mif file for Whole-Brain
%      tractography.
% mask: string, filename for a .mif format of a mask. Use the *_wm.mif file for Whole-Brain
%      tractography.
% mode: Tracking mode: {'prob' | 'stream'} for probabilistic or
%       deterministic tracking. 
% nSeeds: The number of fibers to generate.
% bkgrnd: on unix, whether to perform the operation in another process
% verbose: whether to display standard output to the command window. 
% clobber: Whether or not to overwrite the fiber group if it was already
%          computed
% 
% Franco, Bob & Ariel (c) Vistalab Stanford University 2013 
% Edit GLU 06.2016 added mrTrix versioning

status = 0; results = [];
if notDefined('verbose'),  verbose = false;end
if notDefined('bkgrnd'),    bkgrnd = false;end
if notDefined('clobber'),  clobber = false;end

% Choose the tracking mode (probabilistic or stream)
% mrTrix3 new options
% -algorithm name
%      specify the tractography algorithm to use. Valid choices are: FACT, iFOD1,
%      iFOD2, Nulldist1, Nulldist2, SD_Stream, Seedtest, Tensor_Det, Tensor_Prob
%      (default: iFOD2).
switch mode
  case {'prob','probabilistic tractography'}
    mode_str2 = 'SD_PROB';
    mode_str3 = 'Tensor_Prob';
  case {'stream','deterministic tractogrpahy based on spherical deconvolution'}
    mode_str2 = 'SD_STREAM';
    mode_str3 = 'SD_Stream';
  case {'tensor','deterministic tractogrpahy based on a tensor model'}
    mode_str2 = 'DT_STREAM';
    mode_str3 = 'Tensor_Det';
  otherwise
    error('Input "%s" is not a valid tracking mode', mode); 
end



% In this case there are several changes between mrtrix2 and mrtrix3, I prefer
% to maintain the whole thing for version 2 and create a new one for 3.

if mrtrixVersion == 2
    funcName = 'streamtrack';
    % Generate a UNIX command string.                          
    switch mode_str2
      case {'SD_PROB', 'SD_STREAM'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.csd);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
          strip_ext(mask) , '_', mode, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command
        cmd_str = sprintf('%s %s %s -seed %s -mask %s %s -num %d', ...
                         funcName,mode_str2, files.csd, roi, mask, tck_file, nSeeds); 

      case {'DT_STREAM'}
              % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.dwi);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.dwi), '_' , strip_ext(roi), '_',...
          strip_ext(mask) , '_', mode, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command.
        cmd_str = sprintf('%s %s %s -seed %s -grad %s -mask %s %s -num %d', ...
                 funcName,mode_str2, files.dwi, roi, files.b, mask, tck_file, nSeeds); 
% tckgen iFOD2 
      otherwise
        error('Input "%s" is not a valid tracking mode', mode_str2);
    end
end

if mrtrixVersion == 3
    funcName = 'tckgen';
    % Generate a UNIX command string.                          
    switch mode_str3
      case {'Tensor_Prob', 'SD_Stream'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.csd);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
          strip_ext(mask) , '_', mode, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command
        % Example:
        % tckgen _csd_lmax4.mif        -algo SD_Stream -seed_image _wm.mif               -mask _wm.mif           -num 500000     _csd_lmax4__wm__wm_stream-500000.tck -force
        cmd_str = sprintf('%s %s -algo %s -seed_image %s -mask %s -num %d %s -force', ...
                         funcName, files.csd, ...
                         mode_str3,...
                         roi, ...
                         mask, ...
                         nSeeds, ...
                         tck_file); 

      case {'Tensor_Det'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.dwi);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.dwi), '_' , strip_ext(roi), '_',...
          strip_ext(mask) , '_', mode, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command.
        cmd_str = sprintf('%s %s -algo %s -seed_image %s -grad %s -mask %s -num %d %s -force', ...
                           funcName, files.dwi, ...
                           mode_str3, ...
                           roi, ...
                           files.b, ...
                           mask, ...
                           nSeeds, ...
                           tck_file); 

      otherwise
        error('Input "%s" is not a valid tracking mode', mode_str3);
    end 
end



% Track using the command in the UNIX terminal
if ~(exist(tck_file,'file') ==2)  || clobber == 1
    [status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion);
else
  fprintf('\nFound fiber tract file: %s.\n Loading it rather than retracking',tck_file)
end

% Convert the .tck fibers created by mrtrix to mrDiffusion/Quench format (pdb):
pdb_file = fullfile(pathstr,strcat(strip_ext(tck_file), '.pdb'));
fg = mrtrix_tck2pdb(tck_file, pdb_file);

end 

%%%%%%%%%%%%%
% strip_ext %
%%%%%%%%%%%%%
function [no_ext pathstr] = strip_ext(file_name)
%
% Removes the extension of the files, plus returns the path to the files.
%
[pathstr, no_ext, ext] = fileparts(file_name); 

end