function mrTrixVersion = check_mrTrix_Version
% Check mrTrix version. 
%
% mrTrixVersion = check_mrTrix_Version
%
% mrTrixVersion==0 if mrTrix is not installed. 
% mrTrixVersion==2 or 3 depending on the version
% 
% GLU 06.2016

% Jason, I don't have the AFQ_mrtrix_set_ld_path problem (in linux/OSX), I get rid of it with
% the matlab call itself. Furthermore, it would be circular, since that function
% is using the version information as well. We could us a function call used in
% both versions. Please edit this as it is working for you and I
% will override it in my version if necessary. 


mrTrixVersion = 0;
[status, cmdout] = system('mrconvert -version');
if strfind(cmdout, '== mrconvert 0.2.') % Check if this is how it is written in 0.2
    mrTrixVersion = 2;
end
if strfind(cmdout, '== mrconvert 0.3.')
    mrTrixVersion = 3;
end

