function rgb = AFQ_colormap(cmap,ncolors)

if exist('ncolors','var') && ~isempty(ncolors)
    nsamples = ncolors./4;
else
    nsamples = 64;
    ncolors = 256;
end
if round(nsamples) ~= nsamples
   error('ncolors must be a multiple of 4...sorry') 
end
switch(cmap)
    case{'rgb' 'redblue' 'redgreenblue'}
        rgb = lineargradient([.7 0 0; 1 1 0],nsamples);
        rgb = vertcat(rgb,lineargradient([1 1 0; 0 .7 0],nsamples));
        rgb = vertcat(rgb,lineargradient([0 .7 0; 0 1 1],nsamples));
        rgb = vertcat(rgb,lineargradient([0 1 1; 0 0 .7],nsamples));
    case{'bgr' 'bluered' 'bluegreenred'}
        rgb = lineargradient([0 0 .7; 0 1 1],nsamples);
        rgb = vertcat(rgb,lineargradient([0 1 1; 0 .7 0],nsamples));
        rgb = vertcat(rgb,lineargradient([0 .7 0; 1 1 0],nsamples));
        rgb = vertcat(rgb,lineargradient([1 1 0; .7 0 0],nsamples));
    case{'blueyellow'}
        rgb = lineargradient([0 0 1;1 1 0],256);
end


% figure; colorbar;
% colormap(rgb)
% caxis([1 ncolors])