function afq = AFQ_run_sge(afq, numsubs)
% Run AFQ using a sun grid engine to parallelize the computations
%
% afq = AFQ_run_sge(afq, numsubs)

% set a new outdir
afq.params.outdir = fullfile(afq.params.outdir,'sge');
if ~exist(afq.params.outdir)
    mkdir(afq.params.outdir);
end
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
   % Name the job
   jobname{ii} = sprintf('AFQ%d_%d',ii,round(rand*1000));
   % Process this subject on the grid
   sgerun2('AFQ_run([],[],afq);',jobname{ii},1);
end

% Wait for them to be done
done = 0;
fprintf('\n\n --Waiting for the grid to finish--\n\n')
while done == 0
   % wait 5 seconds
    pause(5) 
    % Check the job status
    [~,stat] = system('qstat');
    % Search for jobs with AFQ in the title
    jobnums = strfind(stat,'AFQ');
    if isempty(jobnums)
        done = 1;
    end
end

%% Reconstruct AFQ structure
% Load up the afq structure for the last subject that was run
load(fullfile(afq.params.outdir,afq.params.outname))
afqfull = afq;
for ii = numsubs
    afqfile = dir(fullfile(afqfull.params.outdir,['afq*_' num2str(ii) '.mat']));
    if ~isempty(afqfile)
    load(fullfile(afq.params.outdir,afqfile.name));
    afqfull = AFQ_CombineSgeRuns(afqfull, afq, ii);
    else
        fprintf('\nNo afq file found for subject %d',ii)
    end
end

% reset the outdir
afqfull.params.outdir = fileparts(afqfull.params.outdir);
% reasign to output
afq = afqfull;
% save
save(fullfile(afq.params.outdir,['AFQ_sge_' date]),'afq');

%% Delete all the jobs to not waste space
for ii = 1:length(jobname)
    delete(fullfile('~/sgeoutput',['job_' jobname{ii} '*']));
end