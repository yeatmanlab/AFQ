%% This tutorial is a work in progress


%% Cleaning a fiber group and visualizing the cleaning proceedure
% See AFQ_removeFiberOutliers
cd /Users/jyeatman/git/AFQ/data/control_01/dti30;
dt = dtiLoadDt6('dt6.mat');
% Load segmented fiber groups
fg = dtiReadFibers('fibers/MoriGroups.mat');
% PUll out the arcuate fasciculus (this is fiber group 19 within the 
% moriGroups,mat)
arc = fg(19);
% Let's take a look at the fibers (only show 50 fibers for times sake)
AFQ_RenderFibers(arc,'numFibers',50);
% Remove fiber outliers
[arc_clean keep]=AFQ_removeFiberOutliers(arc,3,3,30,[],[],5, true);

%% Rendering a fiber group
% See AFQ_RenderFibers
AFQ_RenderFibers(arc_clean, 'numfibers',100);

%% Rebder a fiber group and jitter the shading of each fiber
AFQ_RenderFibers(arc_clean,'numfibers',100,'jittershading',.5)
%% Render a fiber group and Jitter the color of each fiber
AFQ_RenderFibers(arc_clean,'numfibers',100,'jittercolor',.1)

%% Render a fiber group with each point colored based on FA
% See dtiGetValFromFibers
arc_fa = dtiGetValFromFibers(dt.dt6,arc_clean,inv(dt.xformToAcpc),'fa');
rgb = vals2colormap(arc_fa);
AFQ_RenderFibers(arc_clean,'color',rgb,'numfibers',100)

%% 
%% Render a fiber group with each point colored by node number
AFQ
% Cleaned

% Uncleaned


%% Export individual fiber groups as .pdb (Quench format)
