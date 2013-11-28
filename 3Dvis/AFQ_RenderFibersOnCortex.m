function [msh, fdNii, lightH]=AFQ_RenderFibersOnCortex(fg, segmentation, afq, template, fgnums, colormap)

if notDefined('colormap')
    colormap = 'jet';
end
if notDefined('afq')
    afq = [];
end
if notDefined('template')
    template = [];
end
if notDefined('fgnums')
    fgnums = [];
end

%%
if ~isempty(afq) && ~isempty(template)

    defaults = spm_get_defaults;
    params = defaults.normalise.estimate;
    tImg = readFileNifti(template);
        % Allocatge image space
    fdImg = int32(zeros([size(tImg.data) 24]));
    for ii = 1:length(afq.sub_dirs)
        
        dt = dtiLoadDt6(fullfile(afq.sub_dirs{ii},'dt6.mat'));
        fg = AFQ_get(afq,'cleanfibers',ii);
        fgNames = AFQ_get(afq,'fgnames');
        if length(fgNames)>20
            for jj = 21:length(fgNames)
                fg(jj) = dtiReadFibers(AFQ_get(afq,[fgNames{jj} 'path'],ii));
            end
        end
        fg = fg(fgnums);
        if isfield(dt.files,'t1') && exist(dt.files.t1,'file')
            t1=readFileNifti(dt.files.t1);
            
            % Compute normalization
            [sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(t1.data, t1.qto_xyz, template, params);
            % Loop over fibers
            for jj=1:length(fg)
                % Normalize them
                if ~isempty(fg(jj).fibers)
                    fg_sn = dtiXformFiberCoords(fg(jj), invDef);
                    % And add a 1 to each voxel in the group image that
                    % contains a fiber for this subject
                    fdImg(:,:,:,jj) = fdImg(:,:,:,jj)+int32(dtiComputeFiberDensityNoGUI(fg_sn,tImg.qto_xyz,size(tImg.data),1,[],1));
                end
            end
        end
    end
    %% Render the image
    for ii = 1:size(fdImg,4)
        fdImg(:,:,:,ii) = smooth3(fdImg(:,:,:,ii),'gaussian',7);
    end
    % Tack on an extra volume that will mark voxels with no fibers
    fdImg = cat(4,zeros(size(fdImg(:,:,:,1)))+.0000001,fdImg);
    % Find the volume with the highest fiber density in each voxel
    [~,fdMax] = max(fdImg,[],4);
    % Zero out voxels with no fibers
    fdMax = fdMax-1;
    % Make into a nifti volume
    fdNii = segmentation;
    fdNii.data = fdMax;
    fdNii.fname = ['fg_endpoints.nii.gz'];
    clear fdImg
    [p, msh] = AFQ_RenderCorticalSurface(segmentation, 'overlay',fdNii,'boxfilter',3,'thresh',[1 20],'interp','nearest','cmap','wmdevo_colormap');
else
    % Check if the segmentation is binary or is mrVista format
    if length(unique(segmentation.data(:)))>2
        segmentation.data = uint8(segmentation.data==3 | segmentation.data==4);
    end
    for ii = 1:length(fg)
        fdImg(:,:,:,ii) = smooth3(dtiComputeFiberDensityNoGUI(fg(ii),segmentation.qto_xyz,size(segmentation.data),1,[],1),'gaussian',5);
    end
    % Tack on an extra volume that will mark voxels with no fibers
    fdImg = cat(4,zeros(size(fdImg(:,:,:,1)))+.000001,fdImg);
    % Find the volume with the highest fiber density in each voxel
    [~,fdMax] = max(fdImg,[],4);
    % Zero out voxels with no fibers
    fdMax = fdMax-1;
    % Make into a nifti volume
    fdNii = segmentation;
    fdNii.data = fdMax;
    fdNii.fname = ['fg_endpoints.nii.gz'];
    clear fdImg
   [p, msh, lightH] =  AFQ_RenderCorticalSurface(segmentation, 'overlay',fdNii,'boxfilter',3,'thresh',[1 20],'interp','nearest','cmap',colormap);
end