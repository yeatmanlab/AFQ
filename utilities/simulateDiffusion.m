function simulateDiffusion(step,d,myelinLayers,colors,outname)
% Mave a diffusion simulation movie
%
% simulateDiffusion(step,d,myelinLayers,colors,outname)
%
% Inputs:

%% Simulate diffusion
% Tolerance. Distance water can get from axons
tol = .1;
% Set step size in each direction
if notDefined('step')
    step = [.02 .02 .08];
end
Q = eye(3); Q(1,1)=step(1);Q(2,2)=step(2);Q(3,3)=step(3);
% number of iterations of movie at which point to stop rotating
s = [5 45 103];
% trace duration
td = 10;
% number of water molecules
nw = 30;
% layout axons
L =7;
if notDefined('d')
    d=1;
end
[X,Y] = meshgrid(1:d:5,1:d:5);
x = X(:); y = Y(:);
for ii = 1:length(x)
    fg.fibers{ii} = [ones(L,1)*x(ii) ones(L,1)*y(ii) [1:L]'];
end
% Center coordinate of the fibers
m = mean(vertcat(fg.fibers{:}));
% Render the fibers
if notDefined('myelinLayers')
    myelinLayers = 1;
end
% Define colors of axons
if notDefined('colors');
   colors = repmat([.6 0.3 0.2],myelinLayers,1);
end
for ii = 1:myelinLayers
    if ii == 1
        [lh, tubes] = AFQ_RenderFibers(fg,'radius',.2+ii*.03,'color',colors(ii,:),'jittercolor',0,'jittershading',0,'subdivs',20,'camera',[-85 6]);
    else
        % Shorten slightly
        for k = 1:length(fg.fibers)
            fg.fibers{k}(1,3) = fg.fibers{k}(1,3)+.2;
            fg.fibers{k}(end,3) = fg.fibers{k}(end,3)-.2;
        end
        AFQ_RenderFibers(fg,'radius',.3+ii*.03,'color',colors(ii,:),'jittercolor',0,'jittershading',0,'subdivs',30,'camera',[-85 6],'newfig',0);
    end
end
% get all the fiber coordinates
fibcoords = vertcat(fg.fibers{:});
% find the bounding box of the fibers
bb(1,:) = min(fibcoords);
bb(2,:) = max(fibcoords);
set(gca,'color',[0 0 0]);axis('off');set(gcf,'color',[0 0 0]);
AFQ_RenderEllipsoid(Q.*25, [3, -.5, 4], 20,[0 .5 1],0)
axis('image');
camlight(lh,30,-30);
% Get camera view angle
cp = campos;
% Concatenate all fiber coordinates
fc = [];
for ii = 1:length(fg.fibers)
    fc = vertcat(fc,horzcat(tubes.X{ii}(:),tubes.Y{ii}(:),tubes.Z{ii}(:)));
end
% add water randomly distributed aroound center of fibers
C = bsxfun(@plus,randn(nw,3),mean(fc));
% Add as many dimensions as the desired trace length
C = repmat(C,[1 1 td]);
% get figure handle
fh = gcf;
for ii = 1:130
    % Save last iterations positions
    if size(C,3)>1
        for k = size(C,3):-1:2
            C(:,:,k) = C(:,:,k-1);
        end
    end
    % Move the water
    Cnew = computeGaussianDisplacement(Q,squeeze(C(:,:,1)));
    C(:,:,1) = Cnew;
%     % find anything that has left the frame and move it back in
%     del = all(C(:,:,1)<repmat(bb(1,:),size(C,1),1) | C(:,:,1)>repmat(bb(2,:),size(C,1),1)==0,2);
%     C = reshape(C(repmat(del,3,size(C,3))),[sum(del),3,size(C,3)]);
%     np = nearpoints(C(mv==0,:,1)',fibcoords');
%     C(mv==0,:,1) = fibcoords(np,:)+randn(length(np),3).*.05;
    % Loop over each water molecule and render it
    nw = size(C,1);
    for jj = 1:nw
        
        % Render the water at the new position
        wh(jj,1) = AFQ_RenderEllipsoid(eye(3).*.02,C(jj,:,1),10,[0 .5 1],0);
        
        % Plot the new trace
        if size(C,3)>1
            for kk = 2:size(C,3)
                wh(jj,kk) = patchline([C(jj,1,kk);C(jj,1,kk-1)], [C(jj,2,kk);C(jj,2,kk-1)], [C(jj,3,kk);C(jj,3,kk-1)],'edgecolor',[.6 .8 1],'edgealpha',1/(kk-1),'linewidth',3);
            end
        end

    end
    drawnow; % Pause to update frame
    % Add frame to move structure
    mov(ii) = getframe(fh);
    if ii>1 && all(all(mov(ii).cdata(:) == mov(ii-1).cdata(:)))
        keyboard;
    end
    delete(wh(:));
    %campos(cp);
    if ii > 10 && ii < 50
        camorbit(2,2);
        camlight(lh,30,-30);
    elseif ii > 80 && ii < 120
        camorbit(-2,-2);
        camlight(lh,30,-30);
    end
    clear dist ind
end

if notDefined('outname')
    outname = 'diffusionsimulation.avi';
end
movie2avi(mov,outname,'compression','none','fps',8)

return
