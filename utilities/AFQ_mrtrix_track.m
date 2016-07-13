function [status, results, fg, pathstr] = AFQ_mrtrix_track(files, ...
                                                           roi, ...
                                                           mask, ...
                                                           algo, ...
                                                           nSeeds, ...
                                                           bkgrnd, ...
                                                           verbose, ...
                                                           clobber, ...
                                                           mrtrixVersion, ...
                                                           multishell)
%
% function [status, results, fg, pathstr] = mrtrix_track(files, roi, mask, algo, nSeeds, bkgrnd, verbose)
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
% algo: Tracking algorithm: it was 'mode' before. Specify it directly in afq.param.track.mrTrixAlgo 
% nSeeds: The number of fibers to generate.
% bkgrnd: on unix, whether to perform the operation in another process
% verbose: whether to display standard output to the command window. 
% clobber: Whether or not to overwrite the fiber group if it was already
%          computed
% 
% Franco, Bob & Ariel (c) Vistalab Stanford University 2013 
% Edit GLU 06.2016 added mrTrix versioning
% Mode now has been modified with algo, and it is set in the afq structure
% directly in AFQ_Create. Be careful, the options for mrTrix2 and mrTrix3 are
% very different. See AFQ_Create for all the available options.

status = 0; results = [];
if notDefined('verbose'),  verbose = false;end
if notDefined('bkgrnd'),    bkgrnd = false;end
if notDefined('clobber'),  clobber = false;end
if notDefined('mrtrixVersion'),    mrtrixVersion = 3;end
if notDefined('multishell'),  multishell = false;end

% REMOVE THIS IN THE REVIEW
% See AFQ_Create, as indicated by Jason I included the algo in the creation
% process as part of the afq structure. 

% Choose the tracking mode (probabilistic or stream)
% mrTrix3 new options
% -algorithm name
%      specify the tractography algorithm to use. Valid choices are: FACT, iFOD1,
%      iFOD2, Nulldist1, Nulldist2, SD_Stream, Seedtest, Tensor_Det, Tensor_Prob
%      (default: iFOD2).
% switch mode
%   case {'prob','probabilistic tractography'}
%     mode_str2 = 'SD_PROB';
%     mode_str3 = 'Tensor_Prob';
%   case {'stream','deterministic tractogrpahy based on spherical deconvolution'}
%     mode_str2 = 'SD_STREAM';
%     mode_str3 = 'SD_Stream';
%   case {'tensor','deterministic tractogrpahy based on a tensor model'}
%     mode_str2 = 'DT_STREAM';
%     mode_str3 = 'Tensor_Det';
%   otherwise
%     error('Input "%s" is not a valid tracking mode', mode); 
% end



% In this case there are several changes between mrtrix2 and mrtrix3
% Maintaining the whole thing for version 2 and create a new one for 3 so that
% it will be easier in the future to delete the whole mrTrix2 thing

if mrtrixVersion == 2
    funcName = 'streamtrack';
    % Generate a UNIX command string.                          
    switch algo
      case {'SD_PROB', 'SD_STREAM'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.csd);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
          strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command
        cmd_str = sprintf('%s %s %s -seed %s -mask %s %s -num %d', ...
                         funcName,algo, files.csd, roi, mask, tck_file, nSeeds); 

      case {'DT_STREAM'}
              % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.dwi);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.dwi), '_' , strip_ext(roi), '_',...
          strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command.
        cmd_str = sprintf('%s %s %s -seed %s -grad %s -mask %s %s -num %d', ...
                 funcName,algo, files.dwi, roi, files.b, mask, tck_file, nSeeds); 
      otherwise
        error('Input "%s" is not a valid tracking algorithm in mrTrix2', algo);
    end
end

if mrtrixVersion == 3
    % Generate a UNIX command string.                          
    switch algo
      case {'iFOD2'}
        % Build a file name for the tracks that will be generated.
        % THe file name will contain information regarding the files being used to
        % track, mask, csd file etc.
        [~, pathstr] = strip_ext(files.csd);
        tck_file = fullfile(pathstr,strcat(strip_ext(files.csd), '_' , strip_ext(roi), '_',...
          strip_ext(mask) , '_', algo, '-',num2str(nSeeds),'.tck'));

        % Generate the mrtrix-unix command
        % See examples at the end of this file
        if  ~multishell
            cmd_str = ['tckgen ' files.csd ' ' ...
                       '-algo ' algo ' ' ...
                       '-seed_image ' roi ' ' ...
                       '-mask ' mask ' ' ...
                       '-num ' num2str(nSeeds) ' ' ...
                       tck_file ' ' ...
                       '-force'];
        else 
            cmd_str = ['tckgen ' files.csd ' ' ...
                       '-algo ' algo ' ' ...
                       '-seed_image ' roi ' ' ...
                       '-act ' files.tt5 ' ' ...
                       '-num ' num2str(nSeeds) ' ' ...
                       tck_file ' ' ...
                       '-force'];
        end
      case {'FACT', 'iFOD1', 'Nulldist1', 'Nulldist2', 'SD_Stream','Seedtest', 'Tensor_Det', 'Tensor_Prob'}
        error('Wrapper for algo "%s" not implemented yet. Add the case here and remove it from this list. Use the examples below to build your algorithm.', algo);
      otherwise
        error('Input "%s" is not a valid tracking algo in mrTrix3', algo);
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


% tcken options taken from https://github.com/MRtrix3/mrtrix3/blob/master/testing/tests/tckgen
% tckgen dwi.mif   -algo seedtest    -seed_sphere           0,0,4,4  -number 50000 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/seed_sphere.mif 1000
% tckgen dwi.mif   -algo seedtest    -seed_image            mask.mif -number 3888 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/SIFT_phantom_seeds.mif 26
% tckgen dwi.mif   -algo seedtest    -seed_random_per_voxel mask.mif 27 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/SIFT_phantom_seeds.mif 0.5
% tckgen dwi.mif   -algo seedtest    -seed_grid_per_voxel   mask.mif 3 tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif - | testing_diff_data - tckgen/SIFT_phantom_seeds.mif 0.5
% tckgen peaks.mif -algo fact        -seed_rejection        rejection_seed.mif -number 5000 -minlength 4 -mask SIFT_phantom/mask.mif tmp.tck -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif tmp.mif -force && mrstats tmp.mif -mask SIFT_phantom/upper.mif -output mean > tmp1.txt && mrstats tmp.mif -mask SIFT_phantom/lower.mif -output mean > tmp2.txt && testing_diff_matrix tmp1.txt tmp2.txt 30
% tckgen fods.mif  -algo ifod1       -seed_image            mask.mif -act SIFT_phantom/5tt.mif -number 5000 tmp.tck -force && tckmap tmp.tck -template tckgen/act_terminations.mif -ends_only - | mrthreshold - - -abs 0.5 | testing_diff_data - tckgen/act_terminations.mif 0.5
% tckgen dwi.mif   -algo seedtest    -seed_gmwmi            out.mif -act SIFT_phantom/5tt.mif -number 100000 tmp.tck -force && tckmap tmp.tck -template tckgen/gmwmi_seeds.mif - | mrthreshold - - -abs 0.5 | testing_diff_data - tckgen/gmwmi_seeds.mif 0.5
% tckgen peaks.mif -algo fact        -seed_dynamic          fods.mif -mask SIFT_phantom/mask.mif -number 10000 -minlength 4 -initdir 1,0,0 tmp.tck -nthreads 0 -force && tckmap tmp.tck -template SIFT_phantom/dwi.mif tmp.mif -force && mrstats tmp.mif -mask SIFT_phantom/upper.mif -output mean > tmp1.txt && mrstats tmp.mif -mask SIFT_phantom/lower.mif -output mean > tmp2.txt && testing_diff_matrix tmp1.txt tmp2.txt 50
% tckgen fods.mif  -algo ifod2       -seed_image            mask.mif -mask SIFT_phantom/mask.mif -minlength 4 -number 100 tmp.tck -force
% tckgen fods.mif  -algo sd_stream   -seed_image            mask.mif -mask SIFT_phantom/mask.mif -minlength 4 -number 100 -initdir 1,0,0 tmp.tck -force
% tckgen dwi.mif   -algo tensor_prob -seed_image            mask.mif -mask SIFT_phantom/mask.mif -minlength 4 -number 100 tmp.tck -force
% tckgen fods.mif  -algo ifod1       -seed_image            mask.mif -act SIFT_phantom/5tt.mif -backtrack -number 100 tmp.tck -force
% tckgen dwi.mif   -algo tensor_det  -seed_grid_per_voxel   mask.mif 3 -nthread 0 tmp.tck -force && testing_diff_tck tmp.tck tckgen/tensor_det.tck 1e-2
% tckgen dwi.mif   -algo tensor_det  -seed_grid_per_voxel   mask.mif 3 tmp.tck -force && testing_diff_tck tmp.tck tckgen/tensor_det.tck 1e-2
