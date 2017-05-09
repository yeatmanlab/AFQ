function [L_VOF, R_VOF, L_pArc, R_pArc, L_pArc_vot, R_pArc_vot] = ...
                                            AFQ_FindVOF(wholebrainfgPath,...
                                            L_arcuate,...
                                            R_arcuate,...
                                            fsROIdir,...
                                            outdir,...
                                            thresh,...
                                            v_crit, ...
                                            dt, ...
                                            savefiles, ...
                                            arcThresh, ...
                                            parcThresh)
% Segment the VOF from a wholebrain connectome
%
% [L_VOF, R_VOF, L_pArc, R_pArc, L_pArc_vot, R_pArc_vot] = AFQ_FindVOF(wholebrainfgPath,L_arcuate,R_arcuate,fsROIdir,outdir,thresh,v_crit, dt, savefiles, arcThresh, parcThresh)
%
% This function will take in a wholebrain connectome, a segmented arcuate
% fasciculus and a freesurfer segmentation and return the vertical
% occipital fasciculus (VOF).
%
% Inputs:
%
% wholebrainfgPath - A path (or fg structure) for a wholebrain fiber group.
% L_arcuate        - Segmented arcuate fasciculus (left hemisphere). See
%                    AFQ_SegmentFiberGroups
% R_arcuate        - Segmented arcuate fasciculus (right hemisphere).
% fsROIdir         - Path to a directory containing .mat ROIs of each
%                    cortical region that is segmnted by freesurfer. This
%                    means that you must first run freesurfers recon-all on
%                    a t1 weighted image to get a cortical segmentation.
%                    Next use the function:
%                    fs_roisFromAllLabels(fsIn,outDir,type,refT1)
%                    to convert the freesurfer segmentation into ,mat ROIs
% outdir           - This is where all the outputs will be saved
% thresh           - A fiber must travel vertical for a large proportion of
%                    its length. The default values are likely fine
% vcrit            - To find fibers that we can considder vertical, we must
%                    define a threshold of how far a fiber must travel in
%                    the vertical direction. v_crit defines how much
%                    farther a fiber must travel in the vertical direction
%                    compare to other directions (e.g., longitudinal) to be
%                    a candidate for the VOF. The default is 1.3
% dt               - dt6.mat structure or path to the dt6.mat file.
% savefiles        - Logical indicating whether VOF fiber groups should be
%                    saved
%
% Outputs
% L_VOF, R_VOF     - Left and right hemisphere VOF fiber groups
% L_pArc, R_pArc   - Left and right Posterior arcuate fasciculus fiber
%                    groups. The posterior arcuate is another vertical
%                    fiber bundle that marks the anterior extent of the VOF
% L_pArc_vot,      - Some of the posterior arcuate fibers terminate in
% R_pArc_vot         ventral occipitotemporal cortex. We return this subset
%                    of the posterior arcuate as a separate fiber group here.
%
% Copyright Jason D. Yeatman, September 2014. Code released with:
% Yeatman J.D., Weiner K.S., Pestilli F., Rokem A., Mezer A., Wandell B.A.
% (2014). The vertical occipital fasciculus: A forgotten highway. PNAS.

%% Argument checking and parameter setting
% Path to ROIs
if notDefined('fsROIdir')
    fsROIdir = uigetdir([],'Select an ROI directory');
end
% output directory
if notDefined('outdir')
    outdir = uigetdir([],'Select an output directory');
end
% Remove any fiber that doesn't go vertical (positive z) for thresh% of its
% coordinates
if notDefined('thresh')
    thresh = [.95 .6];
end
% Fibers must travel this much farther vertically than other directions
if notDefined('v_crit')
    v_crit = 1.3;
end
% Default is to save fibers
if notDefined('savefiles')
    savefiles = true;
end
% Default is to define VOF as fibers that have fewer than 20 nodes of
% overlap with the arcuate
if notDefined('arcThresh')
    arcThresh = 20;
end
% Default is to consider fibers that are anterior to the posterior arcuate
% as not part of the VOF
if notDefined('parcThresh')
    parcThresh = 1;
end

%% Find vertical fibers
% From the wholebrain fiber group find all the vertical fibers that
% terminate in ventral occipitotemporal cortex (VOT).

[L_fg_vert, R_fg_vert] = AFQ_FindVerticalFibers(wholebrainfgPath,fsROIdir,outdir,thresh,v_crit);

%% Find the posterior arcuate
[L_pArc, R_pArc] = AFQ_Segment_PostArcuate(dt, wholebrainfgPath);

%% Separate VOF from arcuate
L_VOF = dtiNewFiberGroup('L_VOF');
L_pArc_vot = dtiNewFiberGroup('L_posteriorArcuate_vot');
R_VOF = dtiNewFiberGroup('R_VOF');
R_pArc_vot = dtiNewFiberGroup('R_posteriorArcuate_vot');

if ~isempty(L_fg_vert.fibers)
    
    % Make an arcuate fiber density image
    arcFdImg = dtiComputeFiberDensityNoGUI(L_arcuate,dt.xformToAcpc,size(dt.b0));
    % Theshold image at voxels with >n fibers
    arcFdImg = single(arcFdImg>1);
    % Compute the number of nodes that overlap with the arcuate
    fgVals = dtiGetValFromFibers(arcFdImg,L_fg_vert,inv(dt.xformToAcpc));
    fgMvals = cellfun(@(x) sum(x),fgVals);
    % Remove fibers that overlap with the arcuate for more than arcThresh nodes
    L_VOF.fibers = L_fg_vert.fibers(fgMvals<arcThresh);
    L_pArc_vot.fibers =  L_fg_vert.fibers(fgMvals>=arcThresh);
    
    % From the VOF fiber group, remove any fibers that are further anterior
    % than the core of the posterior arcuate.
    if parcThresh == 1
        ymaxL = mean(cellfun(@(x) mean(x(2,:)),L_pArc.fibers));
        L_VOF_keep = cellfun(@(x) all(x(2,:)<ymaxL),L_VOF.fibers);
        % Add fibers that are being removed to pArc_vot
        L_pArc_vot.fibers = vertcat(L_pArc_vot.fibers,L_VOF.fibers(~L_VOF_keep));
        L_VOF.fibers = L_VOF.fibers(L_VOF_keep);
    end
    
    % Clean into a bundle
    L_VOF=AFQ_removeFiberOutliers(L_VOF,4,100,25);
    L_pArc=AFQ_removeFiberOutliers(L_pArc,4,100,25);
else
    L_VOF = [];
    L_pArc_vot = [];
end

%% Repeat for the right hemisphere
if ~isempty(R_fg_vert.fibers)
    
    % Make an arcuate fiber density image
    arcFdImg = dtiComputeFiberDensityNoGUI(R_arcuate,dt.xformToAcpc,size(dt.b0));
    % Theshold image at voxels with >2 fibers
    arcFdImg = single(arcFdImg>2);
    % Compute the number of nodes that overlap with the arcuate
    fgVals = dtiGetValFromFibers(arcFdImg,R_fg_vert,inv(dt.xformToAcpc));
    fgMvals = cellfun(@(x) sum(x),fgVals);
    % Remove fibers that overlap with the arcuate for more than 20 nodes
    R_VOF.fibers = R_fg_vert.fibers(fgMvals<arcThresh);
    R_pArc_vot.fibers =  R_fg_vert.fibers(fgMvals>=arcThresh);
    
    % From the VOF fiber group, remove any fibers that are further anterior
    % than the core of the posterior arcuate.
    if parcThresh == 1
        ymaxR = mean(cellfun(@(x) mean(x(2,:)),R_pArc.fibers));
        R_VOF_keep = cellfun(@(x) all(x(2,:)<ymaxR),R_VOF.fibers);
        % Add fibers that are being removed to pArc_vot
        R_pArc_vot.fibers = vertcat(R_pArc_vot.fibers,R_VOF.fibers(~R_VOF_keep));
        R_VOF.fibers = R_VOF.fibers(R_VOF_keep);
    end
    
    % Clean into a bundle
    R_VOF=AFQ_removeFiberOutliers(R_VOF,4,100,25);
    R_pArc=AFQ_removeFiberOutliers(R_pArc,4,100,25);
else
    R_VOF=[];
    R_pArc_vot=[];
end

%% Save VOF fibers
if savefiles == 1
    sub_dir = fileparts(dt.dataFile);
    dtiWriteFiberGroup(L_VOF,fullfile(sub_dir,'fibers','L_VOF'));
    dtiWriteFiberGroup(R_VOF,fullfile(sub_dir,'fibers','R_VOF'))
end

