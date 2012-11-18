%% Simulate diffusion
% layout axons
L =100;
[X,Y] = meshgrid(1:5,1:5);
x = X(:); y = Y(:);
for ii = 1:length(x)
    fg.fibers{ii} = [ones(L,1)*x(ii) ones(L,1)*y(ii) [1:L]'];
end
AFQ_RenderFibers(fg,'radius',.3,'color',[.6 0.3 0.2],'jittercolor',0,'jittershading',0,'subdivs',10)
set(gca,'color',[0 0 0]);axis('off');set(gcf,'color',[0 0 0]);
% add water
C = rand(100,3)*3;Q = eye(3);
for ii = 1:100
    C = computeGaussianDisplacement(Q,C);
    for jj = 1:100
        AFQ_RenderEllipsoid(eye(3).*.2,C(jj,:),5,[0 .5 1],0);
    end
    pause(.2)
end