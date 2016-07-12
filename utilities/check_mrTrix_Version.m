function mrTrixVersion = check_mrTrix_Version
% Check mrTrix version. 
%
% mrTrixVersion = check_mrTrix_Version
%
% mrTrixVersion==0 if mrTrix is not installed. 
% mrTrixVersion==2 or 3 depending on the version
% 
% GLU 06.2016


mrTrixVersion = 0;
[status2, cmdout2] = system('which streamtrack');
[status3, cmdout3] = system('which tckgen');

if ~status2 && status3
    mrTrixVersion = 2;
elseif status2 && ~status3
    mrTrixVersion = 3;
elseif ~status2 && ~status3
    error('You have both mrTrix2 and mrTrix3 installed, remove one from the path');
else
    mrTrixVersion = 0;
end

