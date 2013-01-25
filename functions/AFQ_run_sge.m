function afq = AFQ_run_sge(afq, numsubs)
% Run AFQ using a sun grid engine to parallelize the computations
%
% afq = AFQ_run_sge(afq, numsubs)

% Do not compute norms
afq = AFQ_set(afq,'computenorms',0);
% Get the number of subjects
if ~exist('numsubs','var') || isempty(numsubs)
    numsubs = 1:AFQ_get(afq,'numsubs');
end

% Loop over the subjects and push the processing to the grid
for ii = numsubs
   % Set which subject to send to the grid
   afq = AFQ_set(afq,'runsubs',ii);
   % Set the name of the saved file
   afq = AFQ_set(afq, 'outname', ['afq_' num2str(ii)]);
   % Process this subject on the grid
   sgerun2('AFQ_run([],[],afq);',sprintf('AFQ%d',ii),1);
end