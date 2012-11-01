function status = check_spm
% Check if spm is in the matlab search path

p = which('spm');
if ~isempty(p)
    status = true;
else
    status = false;
    fprintf('Please download spm and add it to you matlab search path')
end