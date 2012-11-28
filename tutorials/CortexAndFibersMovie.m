%% Rendering videos of the cortex and fiber tracts
%
% This tutorial will show how to render a mesh of the cortical surface and
% underlying fiber tracts and animate them in a video.  You need
% segmentation of the white matter (the higher the resolution the better)
% saved as a binary nifti image and .mat or .pdb files of the fiber tracts.
%
% Copyright Jason D. Yeatman November 2012

%% Load the example data

% Get AFQ directories
[AFQbase AFQdata AFQfunc AFQutil AFQdoc AFQgui] = AFQ_directories;
% The subject's data is saved in the mesh directory
sdir = fullfile(AFQdata,'mesh');
cd(sdir);

% Load up the manual segmentations of the left arcuate and inferior
% longitudinal fasciculus from Quench
arc = dtiLoadFiberGroup('L_Arcuate.pdb');
ilf = dtiLoadFiberGroup('L_ILF.pdb');

% Remove the endpoints of the fibers so they don't show through the cortex.
% This is not necessary but I did not want the fibers to poke through the
% cortex.
for ii = 1:length(arc.fibers)
    arc.fibers{ii} = arc.fibers{ii}(:,5:end-5);
end
for ii = 1:length(ilf.fibers)
    ilf.fibers{ii} = ilf.fibers{ii}(:,5:end-5);
end

%% Render the fibers and cortex

% First we render the first fiber tract in a new figure window
lightH = AFQ_RenderFibers(arc,'color',[0 .5 1],'numfibers',500,'newfig',1);
% Next add the other fiber tract to the same figure window
AFQ_RenderFibers(ilf,'color',[1 .5 0],'numfibers',500,'newfig',0);

% Add a rendering of the cortical surface to the same figure.
% segmentation.nii.gz is a binary image with ones in every voxel that is
% white matter and zeros everywhere else
p = AFQ_RenderCorticalSurface('segmentation.nii.gz', 'boxfilter', 5, 'newfig', 0);

% Delete the light object and put a new light to the right of the camera
delete(lightH);
lightH=camlight('right');
% Turn of the axes
axis('off');

%% Make a video showing the cortex and underlying tracts

% These next lines of code perform the specific rotations that I desired
% and capture each rotation as a frame in the video. After each rotation we
% move the light so that it follows the camera and we set the camera view
% angle so that matlab maintains the same camera distance.
for ii = 1:45
    % Rotate the camera 5 degrees down
    camorbit(0, -2);
    % Set the view angle so that the image stays the same size
    set(gca,'cameraviewangle',8);
    % Move the light to follow the camera
    camlight(lightH,'right');
    % Capture the current figure window as a frame in the video
    mov(ii)=getframe(gcf);
end

% Now rotate along the diagonal axis
for ii = 46:90
    camorbit(2,2);
    set(gca,'cameraviewangle',8);
    camlight(lightH,'right');
    mov(ii)=getframe(gcf);
end

% Rotate to see the top of the brain
for ii = 91:135
    camorbit(0,2);
    set(gca,'cameraviewangle',8);
    camlight(lightH,'right');
    mov(ii)=getframe(gcf);
end

% And rotate back to our starting position
for ii = 136:181
    camorbit(-2,-2);
    set(gca,'cameraviewangle',8);
    camlight(lightH,'right');
    mov(ii)=getframe(gcf);
end
% Adjust the transparence of the cortex so that the underlying fiber tracts
% begin to appear.

a = 1; % Start off with an alpha of 1
for ii = 182:192
    % Slowly make the cortex fade away
    a = a-.03;
    alpha(p,a);
    % Caputure the frame
    mov(ii)=getframe(gcf);
end

% Continue making the cortex slowly fade and rotate the white matter tracts
for ii = 193:373
    a = a -0.015;
    % Don't let the alpha go below 0 (causes an error)
    if a < 0
        a = 0;
    end
    alpha(p,a);
    camorbit(2,0);
    camlight(lightH,'right');
    set(gca,'cameraviewangle',8);
    mov(ii)=getframe(gcf);
end

% Now that we are back to our starting point bring the cortex back
for ii = 374:408
    a = a+.03;
    if a > 1
        a = 1;
    end
    alpha(p,a);
    mov(ii)=getframe(gcf);
end

%% Save the movie as a .avi file to be played by any movie program

movie2avi(mov,'WhiteMatter.avi','compression','none','fps',15)


return

