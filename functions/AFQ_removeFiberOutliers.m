function [fg keep]=AFQ_removeFiberOutliers(fg,maxDist,maxLen,numNodes,M,count,maxIter,show)
% Remove fibers from a fiber group that differ substantially from the
% mean fiber in the group.
%
% [fg keep]=AFQ_removeFiberOutliers(fg,maxDist,maxLen,numNodes,[M = 'mean'],[count = 0],[maxIter = 5], [show = 0])
%
% Inputs:
% fg       = input fiber group structure to be cleaned
% maxDist  = the maximum gaussian distance a fiber can be from the core of
%            the tract and be retained in the fiber group
% maxLen   = The maximum length of a fiber (in stadard deviations from the
%            mean length)
% numNodes = Each fiber will be resampled to have numNodes points
% M        = median or mean to represent the center of the tract
% count    = whether or not to print pruning results to screen
% maxIter  = The maximum number of iterations for the algorithm
% show     = whether or not to show which fibers are being removed in each
%            iteration. If show == 1 then the the fibers that are being 
%            kept and removed will be rendered in 3-D and the user will be 
%            prompted to decide whether continue with the cleaning
%
% Outputs:
% fg       = Cleaned fiber group
% keep     = A 1xN vector where n is the number of fibers in the origional
%            fiber group.  Each entry in keep denotes whether that fiber
%            was kept (1) or removed (0).
%
%  Example:
%  AFQdata = '/home/jyeatman/matlab/svn/vistadata/AFQ/';
%  fg = dtiReadFibers(fullfile(AFQdata,'fibers/L_Arcuate_uncleaned.mat'));
%  maxDist = 4; maxLen = 4;
%  numNodes = 25; M = 'mean'; count = 1; show = 1;
%  [fgclean keep]=AFQ_removeFiberOutliers(fg,maxDist,maxLen,numNodes,M,count,show);
%  % Note that the output fgclean will be equivalent to
%  fgclean = dtiNewFiberGroup; fgclean.fibers = fg.fibers(keep);
%
% Copyright Stanford Vista Team 2011
% Written by Jason D. Yeatman 11/4/2011
%

% whether to compute fiber core with mean or median
if ~exist('M','var') || isempty(M)
    M='mean';
end

% display the number of fibers that were removed in each iteration
if ~exist('count','var') || isempty(count)
    count = 0;
end

% maximum number of iterations
if ~exist('maxIter','var') || isempty(maxIter)
    maxIter = 5;
end

% show the removed fibers or not
if ~exist('show','var') || isempty(show)
    show = 0;
end

% initialize a vector defining the fibers to keep
keep    = ones(length(fg.fibers),1);
keep_prev = keep;
idx0    = find(keep);

% compute the amount of outliers expected given a normal distribution
if exist('normpdf','file')
    maxoutDist  = normpdf(maxDist);   % Requires statistics toolbox
    maxoutLen  = normpdf(maxLen);
else
end

% Create a variable to save the origional fiber group before starting to
% clean it
fg_orig = fg;

% Continue cleaning fibers until there are no outliers
iter=0; cont=1;
while cont ==1 && iter<=maxIter
    % display the number of fibers in the fiber group
    if count==1
        fprintf('\n%s number of fibers: %.0f ',fg.name,length(fg.fibers))
    end
    % keep track of the number of iterations
    iter = iter+1;
    % First remove fibers that are too long
    Lnorm     = AFQ_FiberLengthHist(fg);
    keep1     = Lnorm < maxLen;
    idx1      = find(keep1);
    fg.fibers = fg.fibers(keep1);
    % Represent each node as a 3d gaussian where the covariance matrix is the
    % gaussian distance of each fibers node from the mean fiber
    if length(fg.fibers) > 20
        [SuperFiber, weights] = AFQ_FiberTractGaussian(fg, numNodes, M);
        % Check if any fibers are farther than maxDist from any node.
        keep2 = weights < maxDist;
        keep2 = sum(keep2) == numNodes;
        idx2  = find(keep2);
        % If there are more fibers that are farther than max dist then would be
        % expected in a gaussian distribution then continue iterating
        cont = sum(~keep2) > length(fg.fibers)*maxoutDist;
        % Remove fibers that are more than maxDist from any node
        fg.fibers = fg.fibers(keep2);
        % Keep track of which fibers in the orgional fiber group are being
        % removed
        keep1(idx1) = keep2;
        keep(idx0)  = keep1;
        idx0        = find(keep);
    else
        % If there are more fibers that are farther than max length then would be
        % expected in a gaussian distribution then continue iterating
        cont = sum(~keep1) > length(fg.fibers)*maxoutLen;
    end
    % TODO -- Render the retained and removed fibers for the user to inspect
    if show ==1
        if ~exist('camera','var') || isempty(camera)
            camera = 'sagittal';
        end
        fg_keep = dtiNewFiberGroup('keep','b',[],[],fg_orig.fibers(logical(keep)));
        % fibers that were just removed in this iteration
        new_remove = ~keep & keep_prev;
        fg_remove = dtiNewFiberGroup('keep','b',[],[],fg_orig.fibers(logical(new_remove)));
        % render fibers to keep in blue and to remove in red
        AFQ_RenderFibers(fg_keep,'color',[0 0 1],'camera',camera,'tubes',0);
        AFQ_RenderFibers(fg_remove,'color',[1 0 0],'camera',camera,'newfig',0,'tubes',0);
        % And render the fiber tract core in green
        AFQ_RenderFibers(SuperFiber,'color',[0 1 0],'radius',3,'camera',camera,'newfig',0);
        % ask for users input on whether to keep cleaning
        remove = input('Would you like to remove the fiber outliers shown in red? y/n  ','s');
        % If the user does not like the results of this round of cleaning
        % than we revert back to the fiber group from the previous round
        % and exit out of the cleaning proceedure
        if remove == 0 || strcmpi(remove,'n')
            fg = fg_orig;
            fg.fibers = fg_orig.fibers(logical(keep_prev));
            return
        end
        % get the camera view
        camera = [];
        camera = campos;
        close;
    end
    % Save the keep variable from the current iteration before updating it
    % in the next iteration
    keep_prev = keep;
end
keep = logical(keep);

return