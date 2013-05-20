function renderFiberNodeHeatmap(fg,numNodes,numfibers)
% Render a fiber group
%
%
if ~exist('numfibers','var') || isempty(numfibers)
    numfibers = length(fg.fibers);
end
for ii = 1:length(fg.fibers)
    cmap{ii} = linspace(0,1,length(fg.fibers{ii}))';
end
crgb = vals2colormap(cmap);

% Render the fiber group
AFQ_RenderFibers(fg, 'color',crgb, 'numfibers',numfibers);
