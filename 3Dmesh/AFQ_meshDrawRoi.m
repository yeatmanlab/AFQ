function [coords, indices, bin, msh, spline_points, spline_meshpoints, spline_meshindices] = AFQ_meshDrawRoi(msh, dilate, fh)
%
% [coords, indices, bin, msh, spline_points, spline_meshpoints, spline_meshindices] = AFQ_meshDrawRoi(msh, dilate, fh)
%
%
% Inputs:
%
% msh    - AFQ mesh structure. See AFQ_meshCreate.m
% dilate - How much to dilate each point. Default is 0 which means no
%          dilation
% fh     - figure handle if you already have a mesh window open. Default: open
%          a new figure window
%
% Outputs:
%
% coords            - The coordinates of each point drawn on the mesh
% indices           - Indices for the vertices within msh.tr.vertices
% bin               - A binary vector denoting which vertices in
%                     msh.tr.vertices are within in ROI
% msh               - A mesh is returned with the ROI drawn on it
% spline_points     - Points from a spline fit to the points
% spline_meshpoints - Those spline points projected onto the mesh surface
% spline_meshindices- Indices into msh.tr.vertices that correspond to
%                     spline_meshpoints.
%
% Example:
%
% msh = AFQ_meshCreate;
% [msh.tr.vertices, msh.tr.faces] = read_surf('/mnt/diskArray/projects/freesurfer/fsaverage/surf/lh.white');
% msh.tr.faces = msh.tr.faces + 1;
% msh.tr.FaceVertexCData = repmat([.8 .7 .6],size(msh.tr.vertices,1),1);
% msh.tr = smoothpatch(msh.tr,[],20);
% [coords, indices, bin, msh] = AFQ_meshDrawRoi([],msh, 7)
% AFQ_RenderCorticalSurface(msh);
if ~exist('fh', 'var') || isempty(fh)
    [p,~,lh]=AFQ_RenderCorticalSurface(msh);
    title(sprintf('Press RETURN to exit\nArrows to rotate; < > to zoom'))
    view(180,-90); camlight(lh,'infinite');
    fh = gcf;
    fh = fh.Number;
    ax = gca;
    % get current vertices
    curver = msh.vertex.current;
end
datacursormode on
dcmObj = datacursormode(fh);
set(dcmObj,'SnapToDataVertex','on','Enable','on','DisplayStyle','window');
ii = 0; sp = []; spline_meshpoints = []; currkey= 0;
hold('on');
while currkey==0
    press = waitforbuttonpress;
    
    % Get points off mesh
    if press == 0
        ii = ii+1;
        point = getCursorInfo(dcmObj);
        coords(ii,:) = point.Position
        plot3(coords(ii,1),coords(ii,2),coords(ii,3),'ko','markerfacecolor','k');
        % get index of the coords
        [~,indices(ii)] = ismember(coords(ii,:), msh.tr.vertices, 'rows');
        % Fit a spline to the points
        if ii >1
            cs = cscvn(coords'); spline_points = fnplt(cs)';
            % Map spline points onto mesh
            [spIndices, spbestSqDist] = nearpoints(spline_points',msh.tr.vertices');
            spline_meshpoints = msh.tr.vertices(spIndices,:,:);
            delete(sp);
            sp = plot3(spline_meshpoints(:,1), spline_meshpoints(:,2), spline_meshpoints(:,3),...
                '-b','linewidth',4);
            % get index of the coords
            [~,spline_meshindices(ii)] =   ismember(spline_meshpoints(ii,:), msh.tr.vertices, 'rows');
        end
        
        % Check for button presses and either rotate camera or exit loop
    elseif press == 1
        keypress=get(gcf,'CurrentKey');
        if strcmp(keypress, 'return') || strcmp(keypress, 'escape')
            currkey=1;
        else
            currkey=0;
        end
        % Rotate camera
        % [az, el] = view;
        if strcmp(keypress, 'rightarrow')
            camorbit(-5, 0); camlight(lh,'infinite');
        elseif strcmp(keypress, 'leftarrow')
            camorbit(5, 0); camlight(lh,'infinite');
        elseif strcmp(keypress, 'uparrow')
            camorbit(0, 5); camlight(lh,'infinite');
        elseif strcmp(keypress, 'downarrow')
            camorbit(0, -5); camlight(lh,'infinite');
        elseif strcmp(keypress, 'period')
            camzoom(ax,1.1);
        elseif strcmp(keypress, 'comma')
            camzoom(ax,0.9);
        elseif strcmp(keypress, 's')
            prompt = 'How much to smooth mesh?';
            sm = input(prompt);fprintf('\nSmoothing %d iterations',sm);
            % Smooth the mesh by the desired number of iterations
            msh = AFQ_meshSet(msh,'vertices',sprintf('smooth%d',sm))
            % delete the old rendering and update
            [az, el] = view; [p, lh, fh, ax, dcmObj, sp] = ...
                rerender(msh, az, el, ii, spline_meshpoints);

        end
    end
    
end

msh = AFQ_meshSet(msh,'vertices',curver);
hold('off');

if exist('dilate' ,'var') && ~isempty(dilate) && dilate>0
    for d = 1:dilate
        for ii = 1:length(indices)
            r = msh.tr.faces==indices(ii);
            tmp = msh.tr.faces(find(sum(r,2)>0),:);
            indices = unique(vertcat(indices,tmp(:)));
        end 
    end
end

% close the figure
close(fh)
msh.tr.FaceVertexCData(indices,:) = repmat([1 0 0],length(indices),1);

% Make a logical vector out of the indices
bin = zeros(size(msh.tr.vertices,1),1);
bin(indices) = 1;

return

function [p, lh, fh, ax, dcmObj, sp] = rerender(msh, az, el, ii, spline_meshpoints)
% Rerender mesh
[p,~,lh]=AFQ_RenderCorticalSurface(msh);
title(sprintf('Press RETURN to exit\nArrows to rotate; < > to zoom'))
view(az, el); camlight(lh,'infinite');
fh = gcf;
fh = fh.Number;
ax = gca;
datacursormode on
dcmObj = datacursormode(fh);
set(dcmObj,'SnapToDataVertex','on','Enable','on','DisplayStyle','window');
hold('on');
if ii >1
    sp = plot3(spline_meshpoints(:,1), spline_meshpoints(:,2), spline_meshpoints(:,3),...
        '-b','linewidth',4);
end
