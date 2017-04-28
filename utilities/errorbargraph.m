function [ph, eh] = errorbargraph(barheights, err, group, color, barwidth, colorby)
% Bar plot with error bars and modifiable colors and grouping
%
% h = errorbargraph(barheights, err, group, color)
% 
% Inputs:
% color   = matrix of rgb values. There must be as many colors as there are
%            bars in a group
% Copyright Jason D. Yeatman, December 2012

figure; hold on;

%% argument checking
if size(barheights,1) > 1
    group = repmat(1:size(barheights,2),size(barheights,1),1);
    barheights = barheights(:);
    % NOTE this will be a problem for nonsymetric bars
    err = err(:);
    group = group(:)
end
if ~exist('group','var') || isempty(group)
    group = 1:length(barheights);
end
if ~exist('barwidth','var') ||isempty(barwidth)
    barwidth = .4;
end
if ~exist('colorby','var') ||isempty(bolorby)
    colorby = 'group'
end
groupNums = unique(group);
if size(groupNums,1)>1
    groupNums = groupNums';
end
ngroups = length(groupNums);
barHwidth = barwidth*0.5;
% Make errors a row vector
if size(err,1) > size(err,2)
    err = err';
end
% Reflect errors if there is only 1 per bar
if size(err,1) == 1
    err(2,:) = err;
    % If there is only 1 per point we assume it is error not a confidence
    % interval
    errtype = 'err';
else
    % If there are 2 we assume that it is a confidence interval
    errtype = 'ci'
end
%% plotting
c = 0;
for ii = groupNums
    idx = find(group == ii);
    nsubs = length(idx);
    % xpos = ((1./nsubs)-0.5*(1./nsubs) : 1./nsubs : 1-(0.5*(.5./nsubs))) +ii-0.5
    xpos = ii-((nsubs./2)*barwidth) : barwidth : ii+((nsubs./2)*barwidth) %ii+((nsubs-1)*barwidth)
    for k = 1:nsubs
        c = c+1;
        X = [xpos(k), xpos(k+1), xpos(k+1), xpos(k), xpos(k)]
        Y = [0, 0, barheights(idx(k)), barheights(idx(k)), 0];
        switch(colorby)
            case 'group'
                ph(c)=patch(X,Y,color(k,:));
            case 'bar'
                ph(c)=patch(X,Y,color(c,:));
        end
        % plot the errorbar
        eX = [xpos(k)+barHwidth, xpos(k)+barHwidth];
        switch(errtype)
            case {'error' 'err'}
                eY = [barheights(idx(k))-err(2,idx(k)),barheights(idx(k))+err(1,idx(k))];
                eh = plot(eX,eY,'-k','linewidth',3);
            case {'ci' 'CI' 'confidenceinterval'}
                eh = plot(eX,err(:,idx(k)),'-k','linewidth',3);
        end
    end
end

%% axes
axis tight