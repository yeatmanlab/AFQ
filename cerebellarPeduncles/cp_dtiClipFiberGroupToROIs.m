function fgOut = cp_dtiClipFiberGroupToROIs(fg,ROI1,ROI2, minDist,type)

%% This function clips the cerebellar peduncles to the two ROIs provided and is based on the original "dtiClipFiberGroupToROIs" (vistasoft).

%  Reorders fibers so that the first point intersects with ROI1 and last point intersects with ROI2.
%  Only the two endpoints will connect to the ROIs. If there are multiple path segments that loop between the ROIs the
%  segment closest to the start of the path will be chosen.
 
%  The following change was made to the original "dtiClipFiberGroupToROIs": 
%  1. 'Type' was added as an input variable to define if the ROIs are 'not' or 'and'

%  Original vistasoft function written by Anthony Sherbondy - April 2008
%  Adapted by Maya Yablonski - May 2016 
%  Edited for publication - January 2019

%  Copyright Stanford Vista Team 2008

%% Input variables:
%  fg        - Fiber group (fg) structure to be clipped, here: the cerebellar peduncles.
%  ROI1      - First waypoint ROI used for clipping the fiber group
%              (needs to be a .mat file that is saved in each subject's ROI directory)
%  ROI2      - Second waypoint ROI used for clipping the fiber group
%              (needs to be a .mat file that is saved in each subject's ROI directory)
%  minDist   - A scalar value defining the minimal distance from a fiber to an ROI 
%              in order to count as "a fiber crossing the ROI" (default: .87)
%  type      - NEW: The type of intersection with the two waypoint ROIs
%              Options: 'and' =
%                       'not' = 

%% Argument checking

if ~exist('minDist', 'var')|| isempty(minDist), minDist=.87; end
if ~exist('type', 'var')|| isempty(type), type = {'and'}; end
    
%% Clip fiber group to ROIs
%  Second intersection was changed to 'last' in order to analyse a larger portion of the fiber group, 
%  i.e. from the beginning of the first ROI to the end of the second ROI.
[fgAnd1, foo, keep1, keep1ID] = cp_dtiIntersectFibersWithRoi([], type, minDist, ROI1, fg,'first');
[fgAnd2, foo, keep2, keep2ID] = cp_dtiIntersectFibersWithRoi([], type, minDist, ROI2, fg,'last');

intList = find(keep1&keep2);
fgOut = dtiNewFiberGroup;
fOutCount = 1;
for ii=1:length(intList)
    origI = intList(ii);
    curFiber = fg.fibers{origI};
    if keep1ID(origI) > keep2ID(origI)
        % Reorient the fibers so that end- and start points keep together
        idvec = keep1ID(origI):-1:keep2ID(origI); . 
    else
        idvec = keep1ID(origI):keep2ID(origI);
    end
    % Only store valid length paths
    if( length(idvec) > 1 )
        fgOut.fibers{fOutCount,1} = curFiber(:,idvec);
        fOutCount = fOutCount+1;
    end
end

return;



