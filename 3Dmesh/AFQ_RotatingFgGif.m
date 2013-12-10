function AFQ_RotatingFgGif(fg, colors, outfile, im, slice)
% Make an animated gif of a rotating fiber group

if notDefined('colors')
    colors = jet(length(fg));
end
if notDefined('outfile')
    outfile = fullfile(pwd,'Fibers.gif');
end
% First we render the first fiber tract in a new figure window
lightH = AFQ_RenderFibers(fg(1),'color',colors(1,:),'numfibers',50,'newfig',1);

for ii = 2:length(fg)
    % Next add the other fiber tracts to the same figure window
    AFQ_RenderFibers(fg(ii),'color',colors(ii,:),'numfibers',50,'newfig',0);
end

% Add an image if one was provided
AFQ_AddImageTo3dPlot(im,slice);

% Delete the light object and put a new light to the right of the camera
delete(lightH);
lightH=camlight('right');
% Turn of the axes
axis('off');
axis('image');axis('vis3d');

for ii = 1:180
    % rotate the camera
    camorbit(2,0);
    % Rotate the light with the camera
    lightH = camlight(lightH,'right');
    % Caputure the frame
    mov(ii)=getframe(gcf);
    % Convert the movie frame to an image
    im = frame2im(mov(ii));
    [imind,cm] = rgb2ind(im,256);
    
    % Write out the image into an animated gif
    if ii == 1;
        imwrite(imind,cm,outfile,'gif', 'Loopcount',inf,'delaytime',0);
    else
        imwrite(imind,cm,outfile,'gif','WriteMode','append','delaytime',0);
    end
end