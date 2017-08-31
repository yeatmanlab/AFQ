function afq = AFQ_run_sge(afq, numsubs, v)
% Run AFQ using a sun grid engine to parallelize the computations
%
% afq = AFQ_run_sge(afq, numsubs, v)
%
% Run afq on a compute cluster (or with matlab parallel toolbox). This
% function is brittle and written for specific clusters I have used. Also
% it will not work with 2 groups os subjects (patients and controls), but
% requires only 1 subject group.
%
% Inputs:
% afq     - AFQ structure
% numsubs - Subjects to run. Defaults to all
% v       - Scaler specifying which grid call to use. v=1 for sgerun1,
%           v=2 for sgerun2 (proclus) and v=3 for local cluseter using
%           matlabpool

% which sge command to use
if notDefined('v')
    v = 2;
end

% Do not compute norms
afq = AFQ_set(afq,'computenorms',0);
% Get the number of subjects
if ~exist('numsubs','var') || isempty(numsubs)
    numsubs = 1:AFQ_get(afq,'numsubs');
end

% if an output directory isn't in the afq structure then set one
if isempty(afq.params.outdir)
    afq = AFQ_set(afq,'outdir','~/sge');
    if ~exist(afq.params.outdir,'dir')
        mkdir(afq.params.outdir);
    end
    fprintf('\n No output directory specified. Writing outputs to %s',afq.params.outdir);
    % If outdir was already defined than make an sge directory within
else
    afq.params.outdir = fullfile(afq.params.outdir,'sge');
    if ~exist(afq.params.outdir,'dir')
        mkdir(afq.params.outdir);
    end
end

%% Use sun grid engine (SGE) to run subjects in parallel
% Loop over the subjects and push the processing to the grid
if v == 1 || v == 2
    for ii = numsubs
        % Set which subject to send to the grid
        afq = AFQ_set(afq,'runsubs',ii);
        % Set the name of the saved file
        afq = AFQ_set(afq, 'outname', ['afq_' num2str(ii)]);
        % Name the job
        jobname{ii} = sprintf('AFQ%d_%d',ii,round(rand*1000));
        % Process this subject on the grid
        if v == 2
            sgerun2('AFQ_run([],[],afq);',jobname{ii},1);
        elseif v == 1
            sgerun('AFQ_run([],[],afq);',jobname{ii},1);
        end
    end
    
    % Wait for jobs to be done
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
end

%% Use parallel toolbox to take advantage of local cores
if v == 3
    % If an existing pool is open, then close it
    pool = gcp('nocreate') ;
    delete(pool);
    % Open a new pool
    pool=parpool;
    parfor ii = numsubs
        % Set which subject to send to the grid
        afq_i = AFQ_set(afq,'runsubs',ii);
        % Set the name of the saved file
        afq_i = AFQ_set(afq_i, 'outname', ['afq_' num2str(ii)]);
        % run AFQ
        AFQ_run([],[],afq_i);
    end
    % close the pool
    delete(pool);
end
% This is a hack to deal with differences between sge and matlab parllel
afq  = AFQ_set(afq,'runsubs',numsubs(end))
afq = AFQ_set(afq, 'outname', ['afq_' num2str(numsubs(end))]);

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
% Reset run subs to be the full sample
afqfull = AFQ_set(afqfull,'runsubs',1:AFQ_get(afq,'nsubs'));
% reasign to output
afq = afqfull;
% save
save(fullfile(afq.params.outdir,['AFQ_sge_' date]),'afq');

%% Delete all the jobs to not waste space
if v == 1 || v == 2
    for ii = 1:length(jobname)
        delete(fullfile('~/sgeoutput',['job_' jobname{ii} '*']));
    end
end