function orig_path = AFQ_mrtrix_set_ld_path(mrtrixVersion)
% Set the library path for mrtrix libraries
% Edited by GLU June 2016, added mrtrix3 compatibility

if mrtrixVersion == 2
    funcName = 'csdeconv';
end
if mrtrixVersion == 3
    funcName = 'dwi2fod';
end


% Save origional environment path
orig_path = getenv('LD_LIBRARY_PATH');
% Get mrtrix path
[status, pathstr] = system(['which ' funcName]);
if status~=0
    error('Please install mrtrix')
end
% This will be the path to the mrtrix libraries
pathstr = fullfile(pathstr(1:end-13),'lib');
% set the environment path
setenv('LD_LIBRARY_PATH',pathstr)