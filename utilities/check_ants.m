function [status, ANTSpaths] = check_ants


[stat(1), ANTSpaths{1}] = system('which WarpImageMultiTransform');
[stat(2), ANTSpaths{2}] = system('which antsIntroduction.sh');

% Check if all functions were installed
status = sum(stat) == 0;