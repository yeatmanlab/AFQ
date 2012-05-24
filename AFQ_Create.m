function afq = AFQ_Create(varargin)
% Create AFQ structure and set default parameters
%
%    afq = AFQ_Create(varargin)
%
% Creates a default automated fiber quantification (AFQ) structure.  The
% default fields are put in place and filled with default values.  The
% arguments to the function are in the form 'parameter','value'.  For
% example, if you call the function with
%
%    AFQ_Create('cutoff',[5 95]);  
%
% The afq.params.cutoff parameter will be set to [5 95].
%
% See Also: AFQ_Set AFQ_Get
%
% Examples:
%    afq = AFQ_Create;
%    afq = AFQ_Create('runmode','test');
%
%
% (c) Jason D. Yeatman December 2011
%

%% Attach a structure of values to the afq structure
vals.fa = [];
vals.md = [];
vals.rd = [];
vals.ad = [];
vals.cl = [];
afq.vals = vals;

%% Attach a cell array of subject names to the afq structure
afq.sub_names = {}

%% Attach a vector of subject id numbers to the afq structure
afq.sub_nums = [];

%% Attach a cell array of subject directories to the afq structure
afq.sub_dirs = {};

%% Attach a vector to define each subject's group
afq.sub_group = {};

%% Attach the tract profile structure to the afq structure

% The names of each tract
fgNames = {'L_ATR' 'R_ATR' 'L_CST' 'R_CST' 'L_Cingulum_C' 'R_Cingulum_C'...
'L_Cingulum_H' 'R_Cingulum_H' 'Callosum_Post' 'Callosum_Ant' 'L_IFOF'...
'R_IFOF' 'L_ILF' 'R_ILF' 'L_SLF' 'R_SLF' 'L_Uncinate' 'R_Uncinate'...
'L_Arcuate' 'R_Arcuate'};
TractProfiles = struct;
% Add a field to afq.TractProfiles for each tract defined above
for ii = 1:length(fgNames)
    afq.TractProfiles.(fgNames{ii}) = AFQ_CreateTractProfile('name',fgNames{ii});
end

%% Set the afq.params structure with default parameters
%  cutoff: The percentile cutoff to be used to determine what is "abnormal"
%  The default is cutoff = [10 90] meaning that subjects who fall below the
%  10th percentile or above the 90th percentile will be considered
%  "abnormal" and their data will be plotted with respect to the control
%  population
%  Example: params.cutoff = [10 90];
afq.params.cutoff = [10 90];
%  run_mode: normally leave this blank unless you just want to test this
%  script on a new dataset.  If this is the case then set run_mode = 'test'
%  in which case many fewer fibers will be tracked and segmented to allow
%  for fast testing/debugging of the script
afq.params.run_mode = [];
% cleanFibers = 1 means that once all the fiber groups have been segmented
% they will be cleaned such that any fiber that (1) any fiber that is more
% than params.maxLen standard deviations above the mean fiber length will
% be removed and and (2) any fiber that is more than params.maxDist
% standard deviations from the core of the tract will be removed.  this
% means that fibers groups will be forced to be a compact bundle
afq.params.cleanFibers = 1;
% Remove fibers that are more than maxDist standard deviations from the
% core of the tract
afq.params.maxDist = 4;
% Remove fibers that are more than maxLen standard deviations from the
% mean fiber length
afq.params.maxLen = 3;
% The number of nodes to represent each fiber
afq.params.numberOfNodes = 100;
% Should we analyze the whole length of the fiber group or just the central
% portion spanning between 2 ROIs
afq.params.clip2rois = 0;
% If params.cleanFibers==1, then this will indicate whether to perform the
% cleaning on just the clipped portion of the tract or on the full tract.
% This may be helpful for tracts like the ILF that are extremely messy with
% looping fibers etc.
afq.params.cleanClippedFibers = 0;
% The directory to save all output figures and results
afq.params.outdir = [];
% Save figures? yes or no
afq.params.savefigs = 0;

% TODO:
%  Write a parameter translation routine based on mrvParamFormat()
%  This will translate all of the parameters into the right format for the
%  afq parameters structure.

%% Modify default parameters based on user input
afq.params = mrVarargin(afq.params, varargin);


%%


return
