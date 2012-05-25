function [AFQbase AFQdata AFQfunc AFQutil AFQdoc AFQgui] = AFQ_directories
% Returns a path to the AFQ directories
%
%   [AFQbase AFQdata AFQfunc AFQutil AFQdoc AFQgui] = AFQ_directories
%
% Written by Jason D. Yeatman, May 2012

AFQfunc = fileparts(which('AFQ_run'));
AFQbase = fileparts(AFQfunc);
AFQdata = fullfile(AFQbase,'data');
AFQutil = fullfile(AFQbase,'utilities');
AFQdoc  = fullfile(AFQbase,'doc');
AFQgui  = fullfile(AFQbase,'gui');