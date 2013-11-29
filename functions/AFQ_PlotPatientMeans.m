function AFQ_PlotPatientMeans(afq_patient,afq_controls,valname,nodes,outdir,varargin)
% Plot patient data against controls
%
%
% Example:
%
% load /biac4/wandell/data/WH/analysis/AFQ_WestonHavens_Full.mat
% afq_controls = afq;
% load /biac4/wandell/data/WH/kalanit/PS/AFQ_PS.mat
% afq_patient = afq;
% AFQ_PlotPatientMeans(afq_patient,afq_controls,'T1_map_lsq_2DTI',[],'age', [53 73])
% AFQ_PlotPatientMeans(afq_patient,afq_controls,'fa',[],'age', [53 73])
% AFQ_PlotPatientMeans(afq_patient,afq_controls,'md',[],'age', [53 73])

% Which nodes to analyze
if notDefined('nodes')
    nodes = 20:80;
end

% Define output directory
if notDefined('outdir')
    outdir = pwd;
end

% Define whether plots should be in one window or separate
subplots = 0;

% Open a new figure window for each plot
f1 = figure; hold('on');
f2 = figure; hold('on');
% Get rid of underscores in the valname for the sake of axis labels
valtitle = valname;
sp = strfind(valname,'_');
if ~isempty(sp)
    valtitle(sp) = ' ';
end
% Get number of fiber groups and their names
nfg = AFQ_get(afq_patient,'nfg');
fgNames = AFQ_get(afq_patient,'fgnames');

% Load the fiber group for the patient
fg_p = AFQ_get(afq_patient,'clean fg',1);

% Set the views for each fiber group
fgviews = {'leftsag', 'rightsag', 'leftsag', 'rightsag', ...
    'leftsag', 'rightsag', 'leftsag', 'rightsage', 'axial', 'axial',...
    'leftsag', 'rightsag', 'leftsag', 'rightsag',  'leftsag', 'rightsag'...
    'leftsag', 'rightsag', 'leftsag', 'rightsag'};
% Set the colormap and color range for the renderings
cmap = 'jet'; crange = [-4 4];
% Loop over each fiber group
for ii = 1:nfg
    % Get the values for the patient and compute the mean
    vals_p = AFQ_get(afq_patient,fgNames{ii},valname);
    % Remove nodes that are not going to be analyzed
    vals_p = vals_p(nodes);
    vals_pm = nanmean(vals_p);
    
    % Get the value for each control and compute the mean
    vals_c = AFQ_get(afq_controls,fgNames{ii},valname);
    vals_c = vals_c(:,nodes);
    vals_cm = nanmean(vals_c,2);
    
    % Filter the values based on subject characteristics if they were
    % defined
    for kk = 1:2:length(varargin)
        mdata = AFQ_get(afq_controls,'metadata',varargin{kk});
        mrange = varargin{kk+1};
        use(:,(kk+1)./2) = mdata >= mrange(1) & mdata <= mrange(2);
    end
    % If no metadata was defined use all the controls
    if notDefined('use')
        use = ones(length(afq_controls.sub_dirs),1);
    end
    use = all(use,2);
    % Remove values for unwanted subjects
    vals_c = vals_c(use,:); vals_cm(use,:);
    
    % Compute control group mean and sd
    m = nanmean(vals_cm);
    sd = nanstd(vals_cm);
    % Repeate for the tract profile
    m_tp = nanmean(vals_c);
    sd_tp = nanstd(vals_c);
    
    % Plot control group means and sd
    figure(f1);
    x = [ii-.2 ii+.2 ii+.2 ii-.2 ii-.2];
    y1 = [m-sd m-sd m+sd m+sd m-sd];
    y2 = [m-2*sd m-2*sd m+2*sd m+2*sd m-2*sd];
    y3 = [m-3*sd m-3*sd m+3*sd m+3*sd m-3*sd];
    %fill(x,y3, [.8 .8 .8],'edgecolor',[0 0 0]);
    fill(x,y2, [.6 .6 .6],'edgecolor',[0 0 0]);
    fill(x,y1,[.4 .4 .4] ,'edgecolor',[0 0 0]);
    % Define the color of the point for the fiber group based on its zscore
    tractcol = vals2colormap((vals_pm - m)./sd,cmap,crange);
    % Plot patient
    plot(ii, vals_pm,'ko', 'markerfacecolor',tractcol);
    
    %% Plot tract profile
    if subplots == 1
        % Tract profiles can either be separate or within the same figure
        figure(f2); subplot(10,2,ii); hold('on');
    else
        figure; hold('on');
    end
    y1 = [m_tp-sd_tp fliplr(m_tp+sd_tp)];
    y2 = [m_tp-2*sd_tp fliplr(m_tp+2*sd_tp)];
    y3 = [m_tp-3*sd_tp fliplr(m_tp+3*sd_tp)];
    x = [nodes fliplr(nodes)];
    %fill(x, y3, [.8 .8 .8] ,'edgecolor',[0 0 0]);
    fill(x, y2, [.6 .6 .6] ,'edgecolor',[0 0 0]);
    fill(x, y1, [.4 .4 .4] ,'edgecolor',[0 0 0]);
    plot(nodes,m_tp,'-k','linewidth',3);
    % Plot patient
    plot(nodes, vals_p, '-r', 'linewidth',3);
    title([fgNames{ii} ' - ' valtitle]);
    ylabel(valtitle); xlabel('Location');
    
    % Set the figures size etc
    set(gcf,'units','inches','position',[1 1 8 5],'paperpositionmode','auto')
    % Save an image
    print(gcf,'-djpeg', '-r150', sprintf('tp_%d',ii));
        
    %% Render the tract profile
    
    % Check if the fiber group is empty
    if ~isempty(fg_p(ii).fibers)
        % Get the tract profile
        tp = afq_patient.TractProfiles(ii);
        % zscore the values along the tract profile length
        vals_z = (vals_p - m_tp)./sd_tp;
        % Add the zscores to the tract profile structure
        tp = AFQ_TractProfileSet(tp,'vals','zscore',vals_z);
        
        % Remove the unwanted nodes from the tract profile before rendering
        tp.coords.acpc = tp.coords.acpc(:,nodes);
        
        % Render the tract profile
        AFQ_RenderFibers(fg_p(ii),'tractprofile', tp,...
            'val', 'zscore', 'numfibers', 200, 'cmap', cmap, ...
            'crange', crange, 'camera', fgviews{ii}, 'radius', [.5 5]); 
        axis('off');
        
        % Set the figures size etc
        set(gcf,'units','inches','position',[1 1 8 5],'paperpositionmode','auto')
        % Save an image
        print(gcf,'-djpeg', '-r150', sprintf('fg_%d',ii));
        % Close the figure
        close(gcf);
    end
    
end

% Save out the figure of the tract means
figure(f1)
set(gca,'xtick',1:nfg,'xticklabel',fgNames,'xlim',[0 nfg+1]);
rotateXLabels(gca,45);
ylabel(valtitle);
% Set the figures size etc
set(f1,'units','inches','position',[1 1 15 5],'paperpositionmode','auto')
% Save an image
print(f1,'-djpeg', '-r150', sprintf('means'));
