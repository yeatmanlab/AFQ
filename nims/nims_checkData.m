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
    if path == 0; DWI = 0; MRQ = 0; clear path; return; end
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

MRQ = 0;

% Check for the existence of mrqdata
if isfield(mrq, 'inputdata_seir') && isfield(mrq.inputdata_seir,'name')
    [~, numfiles] = size(mrq.inputdata_seir.name);
    if numfiles == 4
        MRQ = MRQ+1;
    end
end

% % Check for the existence of spgr mrqdata
if isfield(mrq, 'inputdata_spgr') && isfield(mrq.inputdata_spgr,'name')
    [~, numfiles] = size(mrq.inputdata_spgr.name);
    if numfiles >= 2
        MRQ = MRQ+1;
    end
end


%% Return the files

files = {};

if DWI == 1
    files.dw = dw;
end

if MRQ >= 1
    files.spgr = mrq.inputdata_spgr.name;
    files.seir = mrq.inputdata_seir.name;
end

%% Set the exit code and exit

decision = MRQ+DWI;

switch decision
    case 3
        thecode = 111; %  1 = go forward
        fprintf('[%s] - Required data found!\n',mfilename);
        mrq.domrq = 1;
    case {2,1} 
        thecode = 222; %  2 = not yet, try later
        fprintf('[%s] - Only partial data was found!\n',mfilename);
    case 0
        thecode = 000; %  0 = no go
        fprintf('[%s] - Required data not found!\n',mfilename);
end

% If the user wanted to exit, then go for it - close matlab and return the
% exit code to the terminal. By default we don't do this.
if doexit == 1
    fprintf('Exiting - status = %d\n',thecode);
    exit(thecode);
end



return


