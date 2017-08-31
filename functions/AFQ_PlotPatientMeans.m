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
    outdir = {pwd};
end
if ~iscell(outdir)
    out_dir{1} = outdir; outdir = out_dir;
end

% Define whether plots should be in one window or separate
subplots = 1;

% If no value was defined than find all the values that match for controls
% and the patient
if notDefined('valname')
    valname = {'fa' 'md' 'rd' 'ad'};
elseif ischar(valname)
    tmp = valname;clear valname
    valname{1} = tmp;
end
% Get rid of underscores in the valname for the sake of axis labels
for v = 1:length(valname)
    valtitle{v} = valname{v};
    sp = strfind(valname{v},'_');
    if ~isempty(sp)
        valtitle{v}(sp) = ' ';
    end
end
% Get number of fiber groups and their names
nfg = AFQ_get(afq_patient,'nfg');
fgNames = AFQ_get(afq_patient,'fgnames');

% Set the views for each fiber group
fgviews = {'leftsag', 'rightsag', 'leftsag', 'rightsag', ...
    'leftsag', 'rightsag', 'leftsag', 'rightsag', 'axial', 'axial',...
    'leftsag', 'rightsag', 'leftsag', 'rightsag',  'leftsag', 'rightsag'...
    'leftsag', 'rightsag', 'leftsag', 'rightsag'};
% Slices to add to the rendering
slices = [-5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0;...
    0 0 -5; 0 0 -5; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0; -5 0 0; 5 0 0];
% Set the colormap and color range for the renderings
cmap = AFQ_colormap('bgr'); crange = [-4 4];

% Loop over the subjects
for s = 1:length(outdir)
    % Load the fiber group for the patient
    fg_p = AFQ_get(afq_patient,'clean fg',s);
    
    % Load up the b0 image for the patient
    dt_p = dtiLoadDt6(AFQ_get(afq_patient,'dt6path',s));
    b0_p = readFileNifti(dt_p.files.b0);
    
    % Make an output directory for this subject if there isn't one
    if ~exist(outdir{s},'dir')
        mkdir(outdir{s});
    end
    fprintf('\nImages will be saved to %s\n',outdir{s});
    
    % Loop over the different values
    for v = 1:length(valname)
        % Open a new figure window for the mean plot
        f1 = figure; hold('on');
        % Make an output directory if it does not exist
        vout = fullfile(outdir{s},valname{v});
        if ~exist(vout,'dir')
            mkdir(vout);
        end
        % Loop over each fiber group
        for ii = 1:nfg
            % Get the values for the patient and compute the mean
            vals_p = AFQ_get(afq_patient,fgNames{ii},valname{v});
            % Remove nodes that are not going to be analyzed and only
            % compute for subject #s
            vals_p = vals_p(s,nodes);
            vals_pm = nanmean(vals_p);
            
            % Get the value for each control and compute the mean
            vals_c = AFQ_get(afq_controls,fgNames{ii},valname{v});
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
                htp = figure; hold('on'); %subplot('position',[.1 .25 .45 .6]);
            else
                figure; hold('on');
            end
            y1 = [m_tp-sd_tp fliplr(m_tp+sd_tp)];
            y2 = [m_tp-2*sd_tp fliplr(m_tp+2*sd_tp)];
            y3 = [m_tp-3*sd_tp fliplr(m_tp+3*sd_tp)];
            x = [nodes fliplr(nodes)];
            %fill(x, y3, [.8 .8 .8] ,'edgecolor',[0 0 0]);
            po(1) = fill(x, y2, [.6 .6 .6] ,'edgecolor',[0 0 0]);
            po(2) = fill(x, y1, [.4 .4 .4] ,'edgecolor',[0 0 0]);
            po(3) = plot(nodes,m_tp,'-k','linewidth',3);
            % Plot patient
            po(4) = plot(nodes, vals_p, '-r', 'linewidth',3);
            title(fgNames{ii},'fontname','times','fontsize',24);
            ylabel(upper(valtitle{v}),'fontname','times','fontsize',24);
            xlabel('Location','fontname','times','fontsize',24);
            set(gca,'fontname','times','fontsize',18);
            axis('tight');
            
            % If a subplot was desired then grab the image frame
            if subplots == 1
                % Set image properties
                set(gcf,'units','inches','position',[1 1 10 7],...
                    'color',[1 1 1],'paperpositionmode','auto');
                % Change line thickness
                set(po(1:2),'linewidth',2);
                set(po(3:4),'linewidth',8)
                % Convert the figure to an image
                [tpim, tpimcmap] = frame2im(getframe(gcf)); close(gcf)
                % Change to the figure with the tract profile in it
                figure(htp); subplot('position',[.05 .1 .4 .8]);
                imshow(tpim, tpimcmap)
            end
            %% Render the tract profile
            
            % Check if the fiber group is empty
            if ~isempty(fg_p(ii).fibers)
                % Get the tract profile for patient s
                tp = afq_patient.TractProfiles(s,ii);
                % zscore the values along the tract profile length
                vals_z = (vals_p - m_tp)./sd_tp;
                % Add the zscores to the tract profile structure
                tp = AFQ_TractProfileSet(tp,'vals','zscore',vals_z);
                
                % Remove the unwanted nodes from the tract profile before rendering
                tp.coords.acpc = tp.coords.acpc(:,nodes);
                
                % Render the tract profile
                AFQ_RenderFibers(fg_p(ii),'tractprofile', tp, 'color', [.8 .7 .6],...
                    'val', 'zscore', 'numfibers', 200, 'cmap', cmap, ...
                    'crange', crange, 'camera', fgviews{ii}, 'radius', [1 5]);
                axis('off');
                % Add an image in
                AFQ_AddImageTo3dPlot(b0_p,slices(ii,:));
                
                % Add a title
                title('Z Score','fontname','times','fontsize',24);
                % Change the font of the colorbar
                set(colorbar,'fontname','times','fontsize',18)
                if subplots == 1
                    % Set image properties
                    set(gcf,'units','inches','position',[1 1 10 7],...
                        'color',[1 1 1],'paperpositionmode','auto');
                    % Convert the figure to an image
                    [fgim, imcmap] = frame2im(getframe(gcf)); close(gcf)
                    % Change to the figure with the tract profile in it
                    figure(htp); subplot('position',[.5 0 .45 .95]);
                    imshow(fgim, imcmap)
                end
                % Set the figures size etc
                set(gcf,'units','inches','position',[1 1 21 8],...
                    'color',[1 1 1],'paperpositionmode','auto');
                % Save an image
                print(gcf,'-dpng', '-r150', fullfile(vout,sprintf('fg_%d',ii)));
                % Close the figure
                close(gcf);
            end
            
        end
        
        % Save out the figure of the tract means
        figure(f1)
        set(gca,'xtick',1:nfg,'xticklabel',fgNames,'xlim',[0 nfg+1],'fontname','times','fontsize',12);
        rotateXLabels(gca,45);
        ylabel(valtitle{v});
        % Set the figures size etc
        set(f1,'units','inches','position',[1 1 8 4],'paperpositionmode','auto')
        % Save an image twice. Once in the value specific directory and once in
        % the base directory
        print(f1,'-dpng', '-r150', fullfile(vout,sprintf('001_means')));
        print(f1,'-dpng', '-r150', fullfile(outdir{s},sprintf('001_means_%s',valname{v})));
        
        close(f1);
    end
    
    %% Save an animated gif of the rotating fiber group
    AFQ_RotatingFgGif(fg_p,[],fullfile(outdir{s},'000_RotatingFibers.gif'),b0_p,[1 0 0]);
    
end
