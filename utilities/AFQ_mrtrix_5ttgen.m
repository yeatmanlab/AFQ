function [status, results] = AFQ_mrtrix_5ttgen(input_file, ...
                                                  tt5_filename, ...
                                                  bkgrnd, ...
                                                  verbose, ...
                                                  mrtrixVersion,...
                                                  tool) 
%
% Convert a nifti image to a mrtrix .mif image.
% Do this: 5ttgen fsl T1.mif 5tt.mif
% To do: explore more options of the software
% Example
% AFQ_mrtrix_5ttgen('T1.mif', '5tt.mif') 
% 
% Sthg strange is going on with cluster tmp folders
% I create a /tmp folder in $HOME for every subject. 
% The script deletes all folders, but it wil maintain the $HOME/tmp 

if notDefined('bkgrnd'),  bkgrnd  = false;end
if notDefined('verbose'), verbose = true; end
if notDefined('mrtrixVersion'), mrtrixVersion = 3; end
if notDefined('tool'), tool = 'fsl'; end

% There were problems in the cluster with the /tmp directory, and there
% have been reports with the same problems. Just in case use a tmp dir in
% the home folder. Everything will be deleted automatically after use. 
% If you want to maintain the tmp folders for visualization add -nocleanup
tmpDir = '~/tmp';
if ~exist(tmpDir, 'dir'), mkdir(tmpDir), end

% fsl or freesurfer can be selected
if strcmp(tool, 'fsl')
    cmd_str = ['5ttgen fsl ' input_file ' ' tt5_filename ' ' ...
               '-nocrop -tempdir ' tmpDir];
else
    cmd_str = ['5ttgen freesurfer ' ...
               '-lut $FREESURFER_HOME/FreeSurferColorLUT.txt ' ...
               input_file   ' ' ...
               tt5_filename ' ' ...
               '-nocrop -tempdir ' tmpDir];
    
end

[status,results] = AFQ_mrtrix_cmd(cmd_str, ...
                                  bkgrnd, ...
                                  verbose, ...
                                  mrtrixVersion); 
                              
                              
tt5_basename = strsplit(tt5_filename,'.');
seed_gmwmi=strcat(tt5_basename{1},'_gmwmi.mif');
cmd_str   = ['5tt2gmwmi -force ' tt5_filename  ' ' seed_gmwmi ];
[status,results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose, mrtrixVersion)
