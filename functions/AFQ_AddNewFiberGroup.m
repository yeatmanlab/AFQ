function afq = AFQ_AddNewFiberGroup(afq,fgName,roi1Name,roi2Name,cleanFibers,computeVals)
% THIS FUNCTION IS STILL BEING DEVELOPED
% Add new fiber groups from any segmentation proceedure to an afq structure
%
% afq = AFQ_AddNewFiberGroup(afq, fgNaem, roi1Name, roi2Name, [cleanFibers = true],[computeVals = true])
%
% By default AFQ_run will segment the 20 fiber groups defined in the Mori
% white matter atlas and save all the relevant information in the afq
% structure. Aditional groups defined manually or with other segmentation
% algorithms can be added to the afq structure with AFQ_AddNewFiberGroup.
% By default this function will take in (1) an afq strucutre that has already
% been passed through AFQ_run (see AFQ_Create and AFQ_run), (2) the name of
% the new fiber group which should be in each subject's fibers directory,
% and (3) the name of two waypoint ROIs that the fiber group passes through
% which should be in each subjects ROIs directory. The necesary fields in
% the afq structure will be populated, TractProfiles will be calculated and
% added to the structure (See AFQ_CreateTractProfile and
% AFQ_ComputeTractProperties) and the new afq structure will be returned.
% Feel free to take an already segmented fiber group and associate it with
% 2 new ROIs in order to define the edges of the profile differently.
%
% Inputs:
% afq        - An afq structure output from AFQ_run (see also AFQ_Create)
% fgName     - A string containing the name of the fiber group. The fiber
%              group should be in each subject's fibers directory (with the
%              same name). This can be either a .mat file or a .pdb file.
%              If fgName is empty then the fiber group will be defined by
%              intersecting the subject's wholeBrainFG with the two ROIs.
% roi1Name   - A string containing the name of the first waypoint ROI. The
%              ROI should be in each subject's ROIs directory 
% roi2Name   - A string containing the name of the second ROI
% computeVals- Logical defining whether TractProfiles should be calculated
%              for the new fiber group for all of the diffusion (and other) 
%              parameters that are  already contained in the afq structure 
%
% Copyright Jason D. Yeatman November 2012