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

%% Load the AFQ test data
try
	fprintf('\nTesting AFQ_directories.m');
    [AFQbase AFQdata AFQfunc AFQutil AFQdoc AFQgui] = AFQ_directories;
catch ME
    ME.stack
end
testData = load(fullfile(AFQdata,'afq.mat'));

%% Test AFQ_Create.m
sub_dirs = {[AFQdata '/patient_01/dti30'], [AFQdata '/patient_02/dti30']...
    [AFQdata '/patient_03/dti30'], [AFQdata '/control_01/dti30']...
    [AFQdata '/control_02/dti30'], [AFQdata '/control_03/dti30']};
    sub_group = [1, 1, 1, 0, 0, 0]; 
try
    fprintf('\nTesting AFQ_Create.m')
	afq = AFQ_Create('run_mode','test', 'sub_dirs', sub_dirs, 'sub_group', sub_group, 'showfigs', 0); 
catch ME
    ME.stack
end

%% Test AFQ_run.m
try
	fprintf('\nTesting AFQ_run.m')
	[afq patient_data control_data norms abn abnTracts] = AFQ_run(sub_dirs, sub_group, afq);
catch ME
	ME.stack
end

%% Compare values to saved data
vals = AFQ_get(afq,'vals','fa');
valsTest = AFQ_get(testData.afq,'vals','fa');
t = (valsTest - vals) < .001;
eq = find(t ==0);
if isempty(eq)
	fprintf('\nNew FA values match saved FA values');
else
	warning('\nNEW VALUES DO NOT MATCH OLD FA VALUES');
end

