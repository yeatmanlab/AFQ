function [status, mrTrixpaths] = check_mrTrix(version)
% Check to see if mrTrix functions are installed on the system
%
% [status mrTrixpaths] = check_mrTrix
%
% status==1 if mr trix is installed. status ==0 if not installed
% path is a cell array of paths to mrTrix functions

% Edited by GLU 06.2016: mrtrix version compatibility: maintain here a complete
% list of all functions changed.  See mrTrix3 Manual for the complete table of
% changes: MRtrix 0.2 equivalent commands


% Beware, the code is being maintained for both mrTrix2 and mrTrix3. 
% Consider that mrTrix2 is called obsolete by the developers, with no updates 
% http://community.mrtrix.org/t/mrtrix-tutorial-error/141
% Function names change, and there are many new options in mrTrix3.
if notDefined('version'), version = 3; end


[stat(1), mrTrixpaths{1}] = system('which mrview');
[stat(2), mrTrixpaths{2}] = system('which mrconvert');
[stat(3), mrTrixpaths{3}] = system('which mrconvert');
[stat(4), mrTrixpaths{4}] = system('which dwi2tensor');
if version == 2
    [stat(5), mrTrixpaths{5}] = system('which estimate_response'); 
    [stat(6), mrTrixpaths{6}] = system('which csdeconv'); 
    [stat(7), mrTrixpaths{7}] = system('which disp_profile');  
    [stat(8), mrTrixpaths{8}] = system('which tensor2FA'); 
    [stat(9), mrTrixpaths{9}] = system('which mrmult'); 
    [stat(10), mrTrixpaths{10}] = system('which tensor2vector'); 
    [stat(11), mrTrixpaths{11}] = system('which streamtrack'); 

end
if version == 3
    [stat(5), mrTrixpaths{5}] = system('which dwi2response'); 
    [stat(6), mrTrixpaths{6}] = system('which dwi2fod'); 
    [stat(7), mrTrixpaths{7}] = system('which shview');
    [stat(8), mrTrixpaths{8}] = system('which tensor2metric'); % Use -fa output option
    [stat(9), mrTrixpaths{9}] = system('which mrcalc'); % Ej. mrcalc A.mif B.mif -mult out.mif
    [stat(10), mrTrixpaths{10}] = system('which tensor2metric'); % Use -vector output option
    [stat(11), mrTrixpaths{11}] = system('which tckgen');
end

% Check if all functions were installed
status = sum(stat) == 0;