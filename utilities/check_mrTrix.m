function [status, mrTrixpaths] = check_mrTrix
% Check to see if mrTrix functions are installed on the system
%
% [status mrTrixpaths] = check_mrTrix
%
% status==1 if mr trix is installed. status ==0 if not installed
% path is a cell array of paths to mrTrix functions

[stat(1), mrTrixpaths{1}] = system('which mrview');
[stat(2), mrTrixpaths{2}] = system('which mrconvert');
[stat(3), mrTrixpaths{3}] = system('which mrconvert');
[stat(4), mrTrixpaths{4}] = system('which dwi2tensor');
[stat(5), mrTrixpaths{5}] = system('which estimate_response');
[stat(6), mrTrixpaths{6}] = system('which csdeconv');

% Check if all functions were installed
status = sum(stat) == 0;