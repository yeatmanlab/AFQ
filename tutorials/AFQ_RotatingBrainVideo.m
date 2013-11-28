function AFQ_RotatingBrainVideo(fg,msh,colors,outfile)
%% Rendering videos of the cortex and fiber tracts
%
% This tutorial will show how to render a mesh of the cortical surface and
% underlying fiber tracts and animate them in a video.  You need
% segmentation of the white matter (the higher the resolution the better)
% saved as a binary nifti image and .mat or .pdb files of the fiber tracts.
%
% Copyright Jason D. Yeatman November 2012


%% Render the fibers and cortex

% First we render the first fiber tract in a new figure window
lightH = AFQ_RenderFibers(fg(1),'color',colors(1,:),'numfibers',250,'newfig',1);

for ii = 2:length(fg)
    % Next add the other fiber tract to the same figure window
    AFQ_RenderFibers(fg(ii),'color',colors(ii,:),'numfibers',500,'newfig',0);
    
end

% Add a rendering of the cortical surface to the same figure.
% segmentation.nii.gz is a binary image with ones in every voxel that is
% white matter and zeros everywhere else
p = AFQ_RenderCorticalSurface(msh, 'boxfilter', 5, 'newfig', 0);

% Delete the light object and put a new light to the right of the camera
delete(lightH);
lightH=camlight('right');
% Turn of the axes
axis('off');

%% Make a video showing the cortex and underlying tracts

% Adjust the transparence of the cortex so that the underlying fiber tracts
% begin to appear.

a = 1; % Start off with an alpha of 1
for ii = 1:10
    % Slowly make the cortex fade away
    a = a-.03;
    alpha(p,a);
    % Caputure the frame
    mov(ii)=getframe(gcf);
end

% Continue making the cortex slowly fade and rotate the white matter tracts
for ii = 11:191
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
for ii = 192:226
    a = a+.03;
    if a > 1
        a = 1;
    end
    alpha(p,a);
    mov(ii)=getframe(gcf);
end

%% Save the movie as a .avi file to be played by any movie program

movie2avi(mov,outfile,'compression','none','fps',15)


return

