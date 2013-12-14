function rgb = AFQ_colormap(cmap)

switch(cmap)
    case{'rgb' 'redblue' 'redgreenblue'}
        nsamples=64;
        rgb = lineargradient([.7 0 0; 1 1 0],nsamples);
        rgb = vertcat(rgb,lineargradient([1 1 0; 0 .7 0],nsamples));
        rgb = vertcat(rgb,lineargradient([0 .7 0; 0 1 1],nsamples));
        rgb = vertcat(rgb,lineargradient([0 1 1; 0 0 .7],nsamples));
    case{'bgr' 'bluered' 'bluegreenred'}
        nsamples=64;
        rgb = lineargradient([0 0 .7; 0 1 1],nsamples);
        rgb = vertcat(rgb,lineargradient([0 1 1; 0 .7 0],nsamples));
        rgb = vertcat(rgb,lineargradient([0 .7 0; 1 1 0],nsamples));
        rgb = vertcat(rgb,lineargradient([1 1 0; .7 0 0],nsamples));
end


figure; colorbar;
colormap(rgb)
caxis([-3 3])