function [mshRoi, vRoi, msh] = AFQ_meshDrawRoi(msh, dilate, voldata, fill_range, volanat, fh, outname)
%
% [mshRoi, vRoi, msh] = AFQ_meshDrawRoi(msh, dilate, voldata, fill_range, fh, outname)
%
%
% Inputs:
%
% msh       - AFQ mesh structure. See AFQ_meshCreate.m
% dilate    - How much to dilate each point. Default is 0 which means no
%             dilation
% voldata   - Nifti with fMRI data (contrast map or p-values)
% fill_range- Range of values to fill in the volume
% fh        - figure handle if you already have a mesh window open. Default: open
%             a new figure window
% outname   - If defined, will save ROIs with this name and path. if
%             outname is supplied as a blank [] then a file gui will open.
%
% Outputs:
%
% mshRoi.coords            - The coordinates of each point drawn on the mesh. If a
%                            volume fill roi then it will be a nifti image with
%                            the dimensions of the volume.
% mshRoi.indices           - Indices for the vertices within msh.tr.vertices
% mshRoi.bin               - A binary vector denoting which vertices in
%                            msh.tr.vertices are within in ROI
% mshRoi.msh               - A mesh is returned with the ROI drawn on it
% mshRoi.spline_points     - Points from a spline fit to the points
% mshRoi.spline_meshpoints - Those spline points projected onto the mesh surface
% mshRoi.spline_meshindices- Indices into msh.tr.vertices that correspond to
%                            spline_meshpoints.
% vRoi                     - Nifti volume ROI
% msh                      - Return the mesh with the ROI drawn on it
%
% Example 1 Draw line on mesh:
%
% msh = AFQ_meshCreate;
% [msh.tr.vertices, msh.tr.faces] = read_surf('/mnt/diskArray/projects/freesurfer/fsaverage/surf/lh.white');
% msh.tr.faces = msh.tr.faces + 1;
% msh.tr.FaceVertexCData = repmat([.8 .7 .6],size(msh.tr.vertices,1),1);
% msh.tr = smoothpatch(msh.tr,[],20);
% [mshRoi,[], msh] = AFQ_meshDrawRoi([],msh, 7)
% AFQ_RenderCorticalSurface(msh);
%
% % Example 2 - Create fMRI based ROI and save as nifti image:
%
% voldata = 'functionalOverlay-20-Dec-2018.nii.gz';
% t1class = 't1_class.nii.gz';
% im = niftiRead(t1class);
% msh = AFQ_meshCreate(t1class);
% fill_range = [2 inf];
% [~,vRoi] = AFQ_meshDrawRoi(msh, [], voldata, fill_range,[], '/Users/jyeatman/Documents/vwfa_roi');

% Check if ROI should be filled in a volume or not
if ~exist('voldata', 'var') || isempty(voldata)
    drawtype = 'line';
else
    drawtype = 'volfill';
end

% If paths were given then load files
if exist('volanat','var') && ~isempty(volanat) && ischar(volanat)
    volanat = niftiRead(volanat);
end
% Check if we should work with an open figure window
if ~exist('fh', 'var') || isempty(fh)
    if strcmp('volfill',drawtype)
        if ischar(voldata), voldata = niftiRead(voldata); end
        voldata.fname = 'voldata.nii.gz';
        msh = AFQ_meshColor(msh,'overlay',voldata, 'thresh',fill_range);
        % Create a binarized volume
        v = voldata.data > fill_range(1) & voldata.data<fill_range(2);
        vRoi = voldata; vRoi.data = zeros(size(v));
    end
    [p,~,lh]=AFQ_RenderCorticalSurface(msh);
    title(sprintf('Press RETURN to exit\nArrows to rotate; < > to zoom; s to smooth'))
    view(180,-90); camlight(lh,'infinite');
    fh = gcf;
    fh = fh.Number;
    ax = gca;
    % get current vertices
    curver = msh.vertex.current;
end

%% Draw ROI

datacursormode on;
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
        coords(ii,:) = point.Position;
        plot3(coords(ii,1),coords(ii,2),coords(ii,3),'ko','markerfacecolor','k');
        % get index of the coords
        [~,indices(ii)] = ismember(coords(ii,:), msh.tr.vertices, 'rows');
        % Fit a spline to the points
        if ii >1
            cs = cscvn(coords'); spline_points = fnplt(cs)';
            % Map spline points onto mesh
            [spline_meshindices, spbestSqDist] = nearpoints(spline_points',msh.tr.vertices');
            spline_meshpoints = msh.tr.vertices(spline_meshindices,:,:);
        end
        
        switch(drawtype)
            case 'volfill'
                imcoords = ceil(mrAnatXformCoords(voldata.qto_ijk, coords(ii,:)));
                vRoi.data = bwselect3(v,imcoords(2),imcoords(1),imcoords(3)) | vRoi.data>0;
                if sum(vRoi.data(:)) == 0
                    fprintf('\nNO VOLUME DATA AT SELECTED VERTEX.\nTRY ANOTHER SPOT\n');
                else
                    roiMsh = AFQ_meshCreate(vRoi,'color',[0 1 0],'smooth',1); patch(roiMsh.tr); shading('interp'); lighting('gouraud');
                end
                % Display on volume anatomy if it exists
                if exist('volanat', 'var') && ~isempty(volanat)
                    if exist('fhanat','var') && isvalid(fhanat), close(fhanat);end
                    fhanat = figure;
                    % Z plane
                    subplot(1,3,1); hold on;
                    if exist('roiMsh','var'),patch(roiMsh.tr); shading('interp'); lighting('gouraud');end
                    AFQ_AddImageTo3dPlot(volanat,[0, 0, coords(ii,3)],[],[],[],[],'overlay',voldata,fill_range(1));hold on
                    view(0,90);axis image; title(sprintf('Z = %.1f',coords(ii,3)));
                    % Plot cursor
                    plot3(coords(:,1),coords(:,2),coords(:,3),'rx','markersize',7);
                    plot3(coords(end,1),coords(end,2),coords(end,3),'rx','markersize',15);
                    
                    % Y plane
                    subplot(1,3,2); hold on;
                    if exist('roiMsh','var'),patch(roiMsh.tr);shading('interp'); lighting('gouraud');end
                    AFQ_AddImageTo3dPlot(volanat,[0, coords(ii,2), 0],[],[],[],[],'overlay',voldata,fill_range(1));
                    view(0,0);axis image; title(sprintf('Y = %.1f',coords(ii,2)));
                    plot3(coords(:,1),coords(:,2),coords(:,3),'rx','markersize',7);
                    plot3(coords(ii,1),coords(ii,2),coords(ii,3),'rx','markersize',15);

                    % X plane
                    subplot(1,3,3); hold on;
                    if exist('roiMsh','var'),patch(roiMsh.tr); shading('interp'); lighting('gouraud');end
                    AFQ_AddImageTo3dPlot(volanat,[coords(ii,1), 0, 0],[],[],[],[],'overlay',voldata,fill_range(1));
                    view(90,0);axis image; title(sprintf('X = %.1f',coords(ii,1)));
                    plot3(coords(:,1),coords(:,2),coords(:,3),'rx','markersize',7);
                    plot3(coords(end,1),coords(end,2),coords(end,3),'rx','markersize',15);
                    % Turn on 3d rotation and move figure window
                    rotate3d; set(fhanat,'position',get(fh,'position')+[-200 -420 500 0])
                    figure(fh);
                end
            case 'line'
                % Draw on the cortical surface if a line ROI
                if ii >1
                    delete(sp);
                    sp = plot3(spline_meshpoints(:,1), spline_meshpoints(:,2), spline_meshpoints(:,3),...
                        '-b','linewidth',4);
                end
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
            camorbit(0, -5); camlight(lh,'infinite');
        elseif strcmp(keypress, 'downarrow')
            camorbit(0, 5); camlight(lh,'infinite');
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

if exist('vRoi', 'var')
    vRoi.data = uint8(vRoi.data);
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

% If the ROI is just one point then splines will be empty
if size(coords,1) == 1
    spline_points = coords;
    spline_meshpoints = coords;
    spline_meshindices = indices;
end
% Put all the variables into a struct
mshRoi = struct('coords', coords, 'indices', indices, 'bin', bin, ...
 'spline_points',spline_points, 'spline_meshpoints', spline_meshpoints, ...
 'spline_meshindices', spline_meshindices);

% Save if desired
if exist('outname', 'var') && ~isempty(outname)
    vRoi.fname = [outname '.nii.gz'];
    niftiWrite(vRoi);
    save(outname, 'mshRoi');
elseif exist('outname', 'var') && isempty(outname)
    [outname,path] = uiputfile('*.nii.gz','File Selection','volRoi.nii.gz');
    vRoi.fname = fullfile(path, outname);
    niftiWrite(vRoi);
    pp = strfind(outname,'.');
    save(fullfile(path, outname(1:pp(1)-1)), 'mshRoi');
end

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
