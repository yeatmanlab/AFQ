function AFQ_test
% Perform unit testing on all the AFQ functions.
% 
% AFQ_test
%
% This is mainly for developers but could be useful for users as well. This
% funciton will perform unit testing on each function within AFQ to make
% sure it (1) runs without an error and (2) returns the correct values.
% This will allow developers to make sure new code they have written will
% not break or change any of the current functionality. Unit testing is
% essential for continuous integration of code. See wiki pages for
% continuous integration:
% http://en.wikipedia.org/wiki/Continuous_integration
% test driven development
% http://en.wikipedia.org/wiki/Test_driven_development
%
% Copyright Jason D. Yeatman November 2012

% Load the AFQ test data
try
    [AFQbase AFQdata AFQfunc AFQutil AFQdoc AFQgui] = AFQ_directories;
    fprintf('\nTesting AFQ_directories.m');
catch ME
    ME.stack
end
load(fullfile(AFQbase,'test','AFQtestdata'));

% Test AFQ_Create
try
    afq = AFQ_Create
    fprint('\nTesting AFQ_Create.m')
catch ME
    ME.stack
end