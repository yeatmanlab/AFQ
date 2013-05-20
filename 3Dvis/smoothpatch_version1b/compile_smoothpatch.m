function compile_smoothpatch
% Compile smoothpatch mex files
%
% compile_smoothpatch

% get current direcotory
wd = pwd;
% get the AFQ directories
AFQbase = AFQ_directories;
sdir = fullfile(AFQbase,'3Dvis','smoothpatch_version1b');
% Change to the smoothpatch directory
cd(sdir);
% Compile mex files
mex smoothpatch_curvature_double.c -v
mex smoothpatch_inversedistance_double.c -v
mex vertex_neighbours_double.c -v
% Switch back to original directory
cd(wd);