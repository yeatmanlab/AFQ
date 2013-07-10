function h = errorbarvectorplot(data,group,colors)


% number of groups
g = unique(group)
% Ignore zeros
g = g(g~=0);

if notDefined('colors')
    colors = hsv(length(g)).*.8;
end
x0 = 1:size(data,2);
for ii = 1:length(g)
    % Pull out this groups data
    vals = data(group==g(ii),:);
    % number of subjects with non nanmeasurements
    n  = sum(~isnan(vals(:,1)));
    % group mean at each node
    m  = nanmean(vals);
    % standard deviation at each node
    sd = nanstd(vals);
    % standard error of the mean at each node
    se = sd./sqrt(n);
    % plot the mean
    h(1,ii) = plot(m,'-','Color',colors(ii,:),'linewidth',1);
    % plot the confidence interval
    
    fill([x0 fliplr(x0)],[m+se fliplr(m-se)],colors(ii,:),'facealpha',.5,'edgealpha',0)
    
    %     h(3,ii) = plot(m+se,'--','Color',colors(ii,:));
    %     h(2,ii) = plot(m-se,'--','Color',colors(ii,:));
end