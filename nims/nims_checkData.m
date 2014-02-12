function nims_checkData(path)
% 
%  status = nims_checkData(path)
% 
% Given a path, returns a status where 1 indicates that we have all the
% files needed to process a particular data set and 0 means that we don't.
% Perhaps 2 will mean that we have only part of the files needed.
% 
% 
% 
% (C) Stanford University, Vista Lab [2014]
% 


%% Input check
if notDefined('path') || ~exist(path,'dir')
    path = uigetdir(pwd,'Choose path');
    if path == 0; DWI = 0; mrQ = 0; clear path; return; end
end



%% Get a list all the relevant DWI files
dw = getDwiFilesStruct(path);

% Create the initial structure
mrQ.name = tempname;
mrQ.RawDir = path;

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

fprintf('Exiting - status = %d\n',thecode);
% exit(thecode);

return


