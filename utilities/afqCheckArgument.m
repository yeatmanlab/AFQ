function [paramOut valOut] = afqCheckArgument(paramIn, valIn)
% Check to make sure afq parameter and value combinations are formatted
% properly
%
% [paramOut valOut] = afqCheckArgument(param, value)
%

switch(paramIn)
    case({'sub_dirs' 'subdirs' 'subjectdirectories'})
        paramOut = 'sub_dirs';
        if ischar(valIn)
            valOut{1} = valIn;
        elseif ~iscell(valIn)
            error('sub_dirs must be a cell array containing paths');
        else
            valOut = valIn;
            paramOut = paramIn;
        end
    otherwise
        paramOut = paramIn;
        valOut   = valIn;
end
