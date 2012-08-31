function status = check_mrvista;
% Check if mrvista is in the matlab search path

s{1} = which('mrDiffusion');
s{2} = which('mrAnatXformCoords');

if ~isempty(s{1}) && ~isempty(s{2})
    status = true;
else
    status = false;
    fprintf('Please download mrVista and add it to your matlab search path');
end