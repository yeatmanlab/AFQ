function AFQ_plot(varargin)
% This function will plot AFQ results in a variety of ways
%
% AFQ_plot(varargin)
%
% Plotting options are described below:
%
% --Group Plots--
%
% AFQ_plot(afq , 'group')
% or
% AFQ_plot('group1 name', data1, 'group2 name', data2, 'group')
% 
% Tract diffusion profiles can be plotted for multiple groups of subjects.
% For each of the 20 tracts the group mean +/- 1 standard error of the mean
% is plotted.  Each group is on the same plot in a different color. The
% propper input to the function is 'group name' as a string followed
% by data which is the structured array of tract diffusion profiles output from
% AFQ_run.  The string 'group' at the end defines that you would like
% a group plot as oppose to one of the other plotting options.
% 
% Examples:
% [patient_data control_data]=AFQ_run(sub_dirs, sub_group);
% AFQ_plot('patient',patient_data,'control',control_data, 'group');
% Options:
% To plot the data only for a subset of the tracts pass in the 'tracts'
% argument followed by a vector of tract numbers to plot
% AFQ_plot('patient',patient_data,'control',control_data,'tracts',[1,2,6],'group');
%
% --Individual Subject Plots--
% 
% AFQ_plot(afq,'individual')
% or
% AFQ_plot(norms, patient_data, 'ci',[10 90], 'legend',subnames, 'individual')
%
% Tract diffusion profiles will be plotted for each patient with respect to
% the control group norms. The 'ci' argument is optional.  This will define
% the width of the outer confidence interval in percentiles to plot for the
% nomrs.  If this argument is not passed in default values of 10th and 90th
% percentiles will be used. The 'legend' argument is also optional and will
% allow a cell array of subject names to be passed in and these will be
% added to a legend on each plot. 
%
% Examples:
% [patient_data control_data]=AFQ_run(sub_dirs, sub_group);
% subnames = {'Sub 1' 'Sub 2' 'Sub 3' 'Sub 4'};
% AFQ_plot(norms,patient_data,'ci',[10 90], 'legend',subnames,'individual');
%
% Copyright Jason D. Yeatman Stanford Vista Team 2011

%% Argument checks and data formatting

% Check if the first argument is an afq structure. If so we will get the
% data from that.
if isafq(varargin{1})
    afq = varargin{1};
    % Get the proper data from the afq structure depending on the plotting
    % arguments
    if sum(strcmpi('group',varargin))
        data{1} = AFQ_get(afq,'control_data');
        data{2} = AFQ_get(afq,'patient_data');
        gnames = {'control', 'patient'};
    elseif sum(strcmpi('individual',varargin))
        data{1} = AFQ_get(afq,'norms');
        data{2} = AFQ_get(afq,'patient_data');
    elseif sum(strcmpi('colormap',varargin)) > 0
        data{1} = AFQ_get(afq,'control_data');   
    end
end

% Divide up the inputs into data and arguments. This is done by looping
% over the number of input arguments and checking if they are structures or
% strings
argnum = 0; datanum = 0;
for ii = 1 : nargin
    if ischar(varargin{ii})
        argnum        = argnum + 1;
        arg{argnum}   = varargin{ii};
    elseif isstruct(varargin{ii}) && ~isafq(varargin{ii})
        datanum       = datanum + 1;
        data{datanum} = varargin{ii};
    end
end

% check if figures should be saved
s = strcmpi('savefigs', varargin);
if sum(s) > 0
    outdir = varargin{find(s)+1};
    savefigs = 1;
else
    savefigs = 0;
end

% Check if the user specified a subset of tracts to plot otherwise plot all
% of them
t = strcmpi('tracts', varargin);
if sum(t) > 0
    tracts = varargin{find(t)+1};
elseif sum(strcmpi('colormap',arg)) > 0
    % For 'colormap' plots the tract data will be in data{1}
    tracts = 1:length(data{1});
elseif sum(strcmpi('group',arg)) > 0 || sum(strcmpi('individual',arg)) > 0;
    % For 'group' and 'individual' plots the tract data will be in data{2}
    tracts = 1:length(data{2});
end
% tracts should be a row vector
if size(tracts,2) < size(tracts,1)
    tracts = tracts';
end

% Check if there is a defined property to plot
p = strcmpi('property', varargin);
if sum(p) > 0
    property = lower(varargin{find(p)+1});
else
    % default property to plot
    property = 'fa';
end

% Check if there is a defined output directory
p = strcmpi('outdir', varargin);
if sum(p) > 0
    outdir = varargin{find(p)+1};
else
    % default property to plot
    outdir = [];
end

% Check if there are defined subjects to plot otherwise plot all subjects
s = strcmpi('subjects', varargin);
if sum(s) > 0
    % For 'colormap' plots the tract data will be in data{1}
    subjects = varargin{find(s)+1};
elseif sum(strcmpi('colormap',arg)) > 0
    subjects = 1:length(data{1}(1).FA(:,1));
elseif sum(strcmpi('group',arg)) > 0 || sum(strcmpi('individual',arg)) > 0;
    % For 'group' and 'individual' plots the tract data will be in data{2}
    subjects = 1:length(data{2}(1).FA(:,1));
end
% subjects should be a row vector
if size(subjects,2) < size(subjects,1)
    subjects = subjects';
end

% generate random numbers for the figure windowss
fignums = ceil(rand(1,max(tracts)).*10000);

% If an afq structure was passed in then get the fiber group names from it
if exist('afq','var')
    fgNames = AFQ_get(afq,'fgnames');
else
    % Otherwise assume that it is the mori groups and we know their names
    fgNames={'Left Thalmic Radiation','Right Thalmic Radiation','Left Corticospinal','Right Corticospinal', 'Left Cingulum Cingulate', 'Right Cingulum Cingulate'...
        'Left Cingulum Hippocampus','Right Cingulum Hippocampus', 'Callosum Forceps Major', 'Callosum Forceps Minor'...
        'Left IFOF','Right IFOF','Left ILF','Right ILF','Left SLF','Right SLF','Left Uncinate','Right Uncinate','Left Arcuate','Right Arcuate'};
end

%% Group plots
%  Plot the mean +/- 1 std. err. for each data set that is sent in 
if sum(strcmpi('group',arg)) == 1
    % If group names were not pulled from the afq structure then the should
    % have been passed in
    if ~exist('gnames','var')
        gnames = arg(1:length(data));
    end
    % define the colors to be used for each groups plot
    c = lines(length(data));
    % make sure that we have enough fiber group names
    if length(fgNames) < max(tracts)
        fgNames{max(tracts)} = [];
    end
    % Loop over all the tracts
    for jj = 1 : length(tracts)
        % For each tract loop over the number of groups and plot each
        % on the same plot
        for ii = 1 : length(data)
            figure(fignums(jj)); hold on;
            % collect the value of interest
            switch(property)
                case 'fa'
                    vals = data{ii}(tracts(jj)).FA;
                    label = 'Fractional Anisotropy';
                case 'rd'
                    vals = data{ii}(tracts(jj)).RD;
                    label = 'Radial Diffusivity';
                case 'ad'
                    vals = data{ii}(tracts(jj)).AD;
                    label = 'Axial Diffusivity';
                case 'md'
                    vals = data{ii}(tracts(jj)).MD;
                    label = 'Mead Diffusivity';
            end
            % number of subjects with measurements for this tract
            n  = sum(~isnan(vals(:,1)));
            % group mean diffusion profile
            m  = nanmean(vals);
            % standard deviation at each node
            sd = nanstd(vals);
            % standard error of the mean at each node
            se = sd./sqrt(n);
            % plot the mean
            h(ii) = plot(m,'-','Color',c(ii,:),'linewidth',3);
            % plot the confidence interval
            plot(m+se,'--','Color',c(ii,:));
            plot(m-se,'--','Color',c(ii,:));
            % label the axes etc.
            xlabel('Location','fontName','Times','fontSize',12);
            ylabel(label,'fontName','Times','fontSize',12);
            title(fgNames{tracts(jj)},'fontName','Times','fontSize',12);
            set(gca,'fontName','Times','fontSize',12);
        end
        % add a legend to the plot
        legend(h,gnames);
    end    
end
%% Individual plots
%  Plot the norms and confidence intervals then add each individual patient
%  to the plot
if sum(strcmpi('individual',arg)) == 1
    % The first set of data will be norms. Collect the property of interest
    norms = data{1};
    % number of nodes to be plotted
    nnodes = size(data{1}.meanFA,1);
    % collect the value of interest
    switch(property)
        case 'fa'
            Meanvals = norms.meanFA;
            SDvals   = norms.sdFA;
            axisScale = [1 nnodes .2 .9];
            label = 'Fractional Anisotropy';
        case 'rd'
            Meanvals = norms.meanRD;
            SDvals   = norms.sdRD;
            axisScale = [1 nnodes .2 .9];
            label = 'Radial Diffusivity';
        case 'ad'
            Meanvals = norms.meanAD;
            SDvals   = norms.sdAD;
            axisScale = [1 nnodes 1 2.1];
            label = 'Axial Diffusivity';
        case 'md'
            Meanvals = norms.meanMD;
            SDvals   = norms.sdMD;
            axisScale = 'auto';
            label = 'Mead Diffusivity';
        otherwise
            % look for the value in the structure and change the case if
            % needed
            if isfield(norms, ['mean' property]);
                fprintf('\nPlotting %s\n',property);
            elseif isfield(norms, ['mean' upper(property)]);
                property = upper(property);
                fprintf('\nPlotting %s\n',property);
            elseif isfield(norms, ['mean' lower(property)]);
                property = lower(property);
                fprintf('\nPlotting %s\n',property);
            else
                error('Property %s could not be found',property);
            end
            Meanvals = norms.(['mean' property]);
            SDvals = norms.(['sd' property]);
            axisScale = 'auto';
            label = property;
    end
    % make a legend if desired
    if sum(strcmpi('legend',varargin)) > 0
        L = varargin{find(strcmpi('legend',varargin))+1};
    else
        L = [];
    end
    % Define the confidence intervals to be plotted. Percentiles are
    % converted to z scores
    if sum(strcmpi('ci',varargin)) > 0
        cutoff = varargin{find(strcmpi('ci',varargin))+1};
    else
        % if no confidence interval is defined use 10 and 90
        cutoff = [10 90];
    end
    cutZ=norminv(cutoff*.01);
    cutZ2=norminv([.25 .75]);
    % plot the norms for each tract
    for jj = tracts
        figure(fignums(jj));hold on;
        % plot 10 and 90 percentile bands unless the user defines another
        % confidence interval
        x = [1:nnodes fliplr(1:nnodes)];
        y = vertcat(Meanvals(:,jj)+max(cutZ)*SDvals(:,jj), flipud(Meanvals(:,jj)+min(cutZ)*SDvals(:,jj)));
        fill(x,y, [.6 .6 .6]);
        clear y
        % plot the 25 and 75 percentile bands
        y = vertcat(Meanvals(:,jj)+max(cutZ2)*SDvals(:,jj), flipud(Meanvals(:,jj)+min(cutZ2)*SDvals(:,jj)));
        fill(x,y, [.3 .3 .3]);
        clear y
        % plot the mean
        plot(Meanvals(:,jj),'-','Color','k','LineWidth',5);
        % scale the axes and name them
        axis(axisScale);
        xlabel('Location','fontName','Times','fontSize',12);
        ylabel(label,'fontName','Times','fontSize',12)
        title(fgNames{jj},'fontName','Times','fontSize',12)
        set(gcf,'Color','w');
        set(gca,'Color',[.9 .9 .9],'fontName','Times','fontSize',12)
    end
    % the second set of data will be individual subjects
    subData = data{2};
    % define the colors to be used for each individual
    c = hsv(length(subjects)).*0.6;
    % Loop over all the tracts
    for jj = tracts
        figure(fignums(jj));
        % Collect the property of interest for tract jj
        subVals = [];
        switch(property)
            case 'fa'
                subVals = subData(jj).FA;
            case 'rd'
                subVals = subData(jj).RD;
            case 'ad'
                subVals = subData(jj).AD;
            case 'md'
                subVals = subData(jj).MD;
            otherwise
                subVals = subData(jj).(property);
        end
        % For each tract loop over the number of subjects and plot each
        % on the same plot with the norms
        cnum = 0;
        for ii = subjects
            cnum = cnum+1;
            h(ii) = plot(subVals(ii,:)','-','Color',c(cnum,:),'linewidth',2);
        end
        % add a legend to the plot if desired
        if ~isempty(L)
            legend(h(subjects),L);
        end
    end 
    
    % Save the figures if desired
    if ~isempty(outdir)
        cd(outdir);
         for jj = tracts
            figure(fignums(jj));
            fname = [fgNames{jj} property];
            print(gcf, '-depsc',fname)
         end
    end
end

%% Colormap plots
% Plot each individual subject and the group mean.  The Color of the group 
% mean will vary based on the value at that location on the tract profile.
if sum(strcmpi('colormap',arg)) == 1
    % Check if a legend is desired
    if sum(strcmpi('legend',varargin)) > 0
        L = varargin{find(strcmpi('legend',varargin))+1};
    else 
        L = [];
    end
    % Check if the colormap range is defined
    if sum(strcmpi('cmap',varargin)) > 0
        FArange = varargin{find(strcmpi('cmap',varargin))+1};
    else
        FArange = [.3 .6];
    end
    % the first set of data will be individual subjects
    subData = data{1};
    % number of nodes to be plotted
    nnodes = length(subData(1).FA(1,:));
    
    % Loop over all the tracts
    for jj = 1 : length(tracts)
        % collect the property of interest for tract jj
        switch(property)
            case 'fa'
                subVals(:,:,tracts(jj)) = subData(tracts(jj)).FA;
                axisScale = [1 nnodes .2 .9];
                label = 'Fractional Anisotropy';
            case 'rd'
                subVals(:,:,tracts(jj)) = subData(tracts(jj)).RD;
                axisScale = [1 nnodes .2 .9];
                label = 'Radial Diffusivity';
            case 'ad'
                subVals(:,:,tracts(jj)) = subData(tracts(jj)).AD;
                axisScale = [1 nnodes 1 2.1];
                label = 'Axial Diffusivity';
            case 'md'
                subVals(:,:,tracts(jj)) = subData(tracts(jj)).MD;
                axisScale = [1 nnodes .6 1.3];
                label = 'Mean Diffusivity';
            otherwise
                subVals(:,:,tracts(jj)) = AFQ_get(afq, fgNames{tracts(jj)},property);
                label = property;
        end
        
        figure(fignums(jj)); hold on;
        % For each tract loop over the number of subjects and plot each
        % on the same axis
        for ii = 1 : length(subVals(:,1))
            h(ii) = plot(subVals(ii,:,tracts(jj))','-','Color',[.5 .5 .5],'linewidth',1);
        end
        % add a legend to the plot if desired
        if ~isempty(L)
            legend(h(1:length(L)),L);
        end
    end
    
    % To plot the mean we will heatmap each location on the mean
    % profile based on it's FA value and plot a sphere there.
    for jj=1:length(tracts)
        figure(fignums(jj));hold on;
        % compute the mean tract profile
        meanTP = nanmean(subVals(:,:,tracts(jj)));
        % Interpolate the mean profile so the heatmap transitions smoothly
        meanTPinterp = interp1(1:100,meanTP,linspace(1,100,1000));
        % Set the colormap
        c = jet(256);
        % Create a evenly spaced linear mapping between values in the tract
        % profile and colors in the colormap
        lm = linspace(FArange(1), FArange(2),256);
        % Compute the apropriate color for each point on the tract profile
        for k = 1:length(meanTPinterp)
            d = [];
            d = abs(meanTPinterp(k) - lm);
            [tmp, mcolor(k)] = min(d);
        end
        % Plot each point on the tract profile a circle of the apropriate
        % color
        for k = 1:1000
            plot(k./10,meanTPinterp(k),'.','Color',c(mcolor(k),:),'markersize',40);
        end
        % scale the axes and name them
        if exist('axisScale') && ~isempty(axisScale)
            axis(axisScale);
        else
            axis('normal');
        end
        xlabel('Location','fontName','Times','fontSize',12);
        ylabel(label,'fontName','Times','fontSize',12)
        title(fgNames{tracts(jj)},'fontName','Times','fontSize',12)
        set(gcf,'Color','w');
        set(gca,'fontName','Times','fontSize',12)
    end
end
%% Save figures if an output directory is defined
if savefigs==1
    cd(outdir)
    for ii=1:length(fa)
        figure(ii);
        set(gcf,'Color','w','InvertHardCopy','off','PaperPositionMode','auto');
        saveas(gcf,['Figure' num2str(ii)],'png');
    end
end

