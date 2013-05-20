function rgbVals = fiberNodeHeatmap(fg)
% Return a cell array of rgb values to heatmap each fiber from node 1:end
%
% rgbVals = fiberNodeHeatmap(fg)

vals = cellfun(@(x) linspace(0,1,length(x))',fg.fibers,'uniformoutput',0);
rgbVals = vals2colormap(vals);