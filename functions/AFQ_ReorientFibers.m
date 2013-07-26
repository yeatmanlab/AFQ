function [fg, flippedfibers]= AFQ_ReorientFibers(fg, roi1, roi2)
% Flip each fiber so that it passes through roi1 before roi2
%
% fg_reoriented = AFQ_ReorientFibers(fg, roi1, roi2)
%
% When computing Tract Profiles it is important that each subject's fiber
% group has the same start and end points. This function will make sure
% this is true by flipping the start and endpoints of any fiber that passes
% through roi2 before roi1. roi1 and roi2 can either be at the endpoints of
% the fiber group or can be waypoint rois along the trajectory of the fiber
% group.
%
% Input:
%
% fg            - Fiber group
% roi1          - The first roi that the fibers should pass through
% roi2          - The second roi that the fibers should pass through
%
% Output:
% fg            - Fiber group with the fibers reoriented appropriately
% flippedfibers - A logical denoting which fibers were flipped
%
% Written by Jason Yeatman and Tingting Liu July 2013

% preallocate a variable to record which fibers are flipped
flippedfibers = false(1,length(fg.fibers));

% Loop over the fibers in a fiber group
for ii = 1:length(fg.fibers)
   % Compute the distance between each fiber node and roi1 and roi2
   [indices_roi1, dist1]= nearpoints(roi1.coords', fg.fibers{ii});
   [indices_roi2, dist2] = nearpoints(roi2.coords', fg.fibers{ii});
   % For this fiber find the index of the point that is closest to roi1 and roi2
   [~,ind1] = min(dist1);
   [~,ind2] = min(dist2);
   % Find the node number that is closest to each roi by using the index
   nearpoint_roi1 = indices_roi1(ind1);
   nearpoint_roi2 = indices_roi2(ind2);
   % We want to flip any fiber which does not pass through roi1 before
   % roi2. Hence we will flip any fiber for which nearpoint_roi1 is not
   % less than nearpoint_roi2
   if nearpoint_roi1 > nearpoint_roi2
       fg.fibers{ii} = fliplr(fg.fibers{ii});
       % Record whether fiber ii was flipped. The default is false.
       flippedfibers(ii) = true;
   end
end

return