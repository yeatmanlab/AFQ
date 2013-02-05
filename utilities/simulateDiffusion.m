%% Simulate diffusion
% Tolerance. Distance water can get from axons
tol = .1;
% Set step size
step = [.02 .02 .06];

% layout axons
L =8;
[X,Y] = meshgrid(1:5,1:5);
x = X(:); y = Y(:);
for ii = 1:length(x)
    fg.fibers{ii} = [ones(L,1)*x(ii) ones(L,1)*y(ii) [1:L]'];
end
% Center coordinate of the fibers
m = mean(vertcat(fg.fibers{:}));
% Render the fibers
[lh, tubes] = AFQ_RenderFibers(fg,'radius',.3,'color',[.6 0.3 0.2],'jittercolor',0,'jittershading',0,'subdivs',10);
set(gca,'color',[0 0 0]);axis('off');set(gcf,'color',[0 0 0]);
% Get camera view angle
cva = camva;
% Concatenate all fiber coordinates
fc = [];
for ii = 1:length(fg.fibers)
    fc = vertcat(fc,horzcat(tubes.X{ii}(:),tubes.Y{ii}(:),tubes.Z{ii}(:)));
end
% add water randomly distributed aroound center of fibers
nw = 100;
C = bsxfun(@plus,rand(nw,3).*3,[1.5 1.5 1.5]);
Q = eye(3); Q(1,1)=step(1);Q(2,2)=step(2);Q(3,3)=step(3);
% get figure handle
fh = gcf;
for ii = 1:100
    % Save last iterations positions
    tmp = C;
    % Move the water
    C = computeGaussianDisplacement(Q,C);
    % Check which water molecules hit fibers
    [~, dist] = nearpoints(C',fc');
    ind = sqrt(dist)<tol;
    % Replace water that is hitting an axon with it's old location
    C(ind,:) = tmp(ind,:);
    for jj = 1:nw
        wh(jj)=AFQ_RenderEllipsoid(eye(3).*.03,C(jj,:),5,[0 .5 1],0);   
    end
    mov(ii) = getframe(fh);
    pause(.2)
    delete(wh);
    camva(cva);
    if ii < 30
        camorbit(3,3);
        camlight(lh,'right');
    end
end