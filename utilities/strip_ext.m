function [no_ext pathstr] = strip_ext(file_name)
%
% Removes the extension of the files, plus returns the path to the files.
%
[pathstr, no_ext, ext] = fileparts(file_name); 

end