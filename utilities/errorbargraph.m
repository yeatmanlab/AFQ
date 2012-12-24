function h = errorbargraph(barheights, err, group, color, barwidth, ax)
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
    err = err(:);
    group = group(:)
end
if ~exist('group','var') || isempty(group)
    group = 1:length(barheights);
end
if ~exist('ax', 'var') || isempty(ax)
    ylim = [min(barheights)-2*max(err) max(barheights)+2*max(err)];
end
groupNums = unique(group)';
ngroups = length(groupNums);
barHwidth = barwidth*0.5;
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
        h(c)=patch(X,Y,color(k,:));
        % plot the errorbar
        eX = [xpos(k)+barHwidth, xpos(k)+barHwidth];
        eY = [barheights(idx(k))-err(idx(k)),barheights(idx(k))+err(idx(k))];
        plot(eX,eY,'-k','linewidth',4);
    end
end

%% axes
set(gca,'ylim',ylim,'xlim',[0 ngroups+1]);