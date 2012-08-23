function status = check_mrTrix
% Check to see if mrTrix functions are installed on the system
%
% status = check_mrTrix
%
% status==1 if mr trix is installed. status ==0 if not installed
%

stat(1) = system('which mrview');
stat(2) = system('which mrconvert');
stat(3) = system('which mrconvert');
stat(4) = system('which dwi2tensor');
stat(5) = system('which estimate_response');
stat(6) = system('which csddeconv');

% Check if all functions were installed
status = sum(stat) == 0