function mrTrixVersion = check_mrTrix_Version
% Check mrTrix version. 
%
% mrTrixVersion = check_mrTrix_Version
%
% mrTrixVersion==0 if mr trix is not installed. 
% mrTrixVersion==2 or 3 depending on the version
% 
% GLU 06.2016

mrTrixVersion = 0;
[status, cmdout] = system('mrconvert -version');
if strfind(cmdout, '== mrconvert 0.2.') % Check if this is how it is written in 0.2
    mrTrixVersion = 2;
end
if strfind(cmdout, '== mrconvert 0.3.')
    mrTrixVersion = 3;
end

