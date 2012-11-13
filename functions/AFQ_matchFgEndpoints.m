function [fg flipped] = AFQ_matchFgEndpoints(fg,startpoint,endpoint)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Reorder each fiber in a group so that they all go the same direction (eg
% anterior-posteror). Each fiber within a fiber group will be reoriented so
% that it's first coordinate is closer to the defined startpoint than it's
% last coordinate is (and vice versa for endpoint).
%
% [fg flipped] = AFQ_matchFgEndpoints(fg,startpoint,endpoint)
%
% Imputs:
% fg        - Fiber group structure containing 1 or more fiber groups
% startpoint- A Nx3 matrix containing the desired startpoint for each fiber
%             in a fiber group. N is the number of fiber groups in fg
% endpoint  - A Nx3 matrix containing the desired endpoints
%
% Outputs:
% fg     - Fiber group with all the fibers flipped
% flipped- Which fibers were flipped
%
% Copyright Jason D. Yeatman, October 2012

%% Argument checking
if ~exist('fg','var') || isempty(fg)
    error('Please supply a fiber group')
end
if ~exist('startpoint','var') || isempty(startpoint)
    startpoint =[];
elseif size(startpoint,2) ~= 3
    startpoint = startpoint';
end
if ~exist('end','var') || isempty(endpoint)
    endpoint=[];
elseif size(endpoint,2) ~=3
    endpoint = endpoint';
end

%% Reorient the fibers within each suplied fiber group
% Number of fiber groups
nFG = length(fg);
for jj = 1:nFG
    % Check if the fiber group is empty
    if isempty(fg(jj).fibers)
        continue
    end
    % Make a function to compute the distance between fiber endpoints and
    % the desired start and endpoint
    distfun = @(x) nearpoints([x(:,1) x(:,end)],[startpoint(jj,:)' endpoint(jj,:)']);
    % Compute distance for each fiber
    [ind, dist] = cellfun(distfun,fg(jj).fibers,'UniformOutput',false);
    % record which fibers need to be flipped
    flipped{jj} = cellfun(@(x) x(1)>x(2), ind,'UniformOutput',false);
    % Flip fibers that need to be flipped CHECK THIS!!!
    fg(jj).fibers = cellfun(@flipfun,fg(jj).fibers,flipped{jj},'UniformOutput',false);
end

return

function fibOut = flipfun(fibIn,ind)
% Function to flip fibers
if(ind ==1)
    fibOut = fliplr(fibIn);
else
    fibOut = fibIn;
end