function [thecode, files] = nims_checkData(path,doexit)
% 
%  status = nims_checkData(path)
% 
% Given a path, returns a status where 111 indicates that we have all the
% files needed to process a particular data set and 000 means we don't.
% 2222 will mean that we have only part of the files needed.
% 
% 
% 
% 
% (C) Stanford University, Vista Lab [2014] - LMPERRY
% 


%% Input check
if notDefined('path') || ~exist(path,'dir')
    path = uigetdir(pwd,'Choose path');
    if path == 0; DWI = 0; mrQ = 0; clear path; return; end
end

if notDefined('exit')
    doexit = 0;
end


%% Get a list all the relevant DWI files
dw = getDwiFilesStruct(path);

% Check for the existence of dwData
if isempty(dw)
    DWI = 0; 
else
    if ( isfield(dw{1},'nifti') && exist(dw{1}.nifti,'file') && ...
            isfield(dw{1},'bval')  && exist(dw{1}.bval,'file')  && ...
            isfield(dw{1},'bvec')  && exist(dw{1}.bvec,'file') )
        DWI = 1;
    else
        DWI = 0;   
    end
end


%% Given mrQ.RawDir, the following function will return all the necessarry
% paths to required nifti files in the mrQ.inputData_spgr and
% mrQ.inputData_seir structures (maybe this should be part of mrQ create).
% Create the initial structure

mrQ.name = tempname;
mrQ.RawDir = path;
mrq = mrQ_initInputData(mrQ);


% Check for the existence of mrqdata
if isfield(mrq.inputdata_seir','name')
    [~, numfiles] = size(mrq.inputdata_seir.name);
    if numfiles == 4
        mrQ = 1;
    end
else
    mrQ = 0;
end


% % Check for the existence of spgr mrqdata
if isfield(mrq.inputdata_spgr','name')
    [~, numfiles] = size(mrq.inputdata_spgr.name);
    if numfiles >= 2
        mrq.domrq = 1;
    end
else
    mrQ = 0;
end


%% Return the files

files = {};
files.dw = dw;
files.spgr = mrq.inputdata_spgr.name;
files.seir = mrq.inputdata_seir.name;


%% Set the exit code and exit

decision = mrQ+DWI;

switch decision
    case 2
        thecode = 111; %  1 = go forward
    case 1
        thecode = 222; %  2 = not yet, try later
    case 0
        thecode = 000; %  0 = no go
end

% If the user wanted to exit, then go for it - close matlab and return the
% exit code to the terminal. By default we don't do this.
if doexit == 1
    fprintf('Exiting - status = %d\n',thecode);
    exit(thecode);
end



return


