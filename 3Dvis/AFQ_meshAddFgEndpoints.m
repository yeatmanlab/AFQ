function msh = AFQ_meshAddFgEndpoints(msh, fg, colors, crange, alpha, weights, distfun, dilate)
% Color mesh vertices based on fiber endpoints
%
% msh = AFQ_meshAddFgEndpoints(msh, fg, colors, crange, alpha, weights, dilate)
% 
% Inputs:
%
% msh
% fg      - Fiber group structure or path to one
% colors  - nx3 matrix of rgb values for which to color fiber endpoint
%           density on the cortical surface.
% crange  - The endpoint density to map to the minumum and maximum color
%           values
% weights - A vector of weights associated with each fiber.
% dilate  - How far along the cortical surface should we dilate the fiber
%           endpoints.
%
% Example:
%
%
% Copyright - Jason D. Yeatman, December 2013

% By default just map each fiber to the nearest mesh vertex
if notDefined('distfun')
    distfun = 'nearpoints';
end

% If weights aren't defined then use a weight of 1 for each fiber
if notDefined('weighted')
    weights = ones(length(fg.fibers),1);
end

% Default to fully opaque coloring
if notDefined('alpha')
    alpha = 1;
end

% jet is default color map
if notDefined('colors')
    colors = jet(256);
end

% Read in the fiber group if a path was given
if ischar(fg)
    fg = fgRead(fg)
end

% By default dilate the endpoints by 1 face along the mesh
if notDefined('dilate')
    dilate = 1;
end

% Get the fiber group name with no spaces and remove any characters that 
% are not allowed for field names
rmchar = {' ','_','1','2','3','4','5','6','7','8','9','0'};
fgname = fg.name;
for ch = 1:length(rmchar)
    if fgname(1)==rmchar{ch}
        fgname = horzcat('x',fgname);
        continue
    end
end
fgname(strfind(fgname,' ')) = '_';
fgname(strfind(fgname,'-')) = '_';

% Get the start and endpoints of the fibers
for ii = 1:length(fg.fibers)
    sp(:,ii) = fg.fibers{ii}(:,1);
    ep(:,ii) = fg.fibers{ii}(:,end);
end

%% Map fiber endpoints to weights on the cortical surface
% There are multiple distance functions that can be used
switch(distfun)
    
    case{'nearpoints'}
        
        % Find the closest mesh vertex to each start coordinate
        msh_indices1 = nearpoints(sp, msh.vertex.origin');
        % Find the closes mesh vertex to each end coordinate
        msh_indices2 = nearpoints(ep, msh.vertex.origin');
        
        % Find the weighted endpoint count at each mesh index. To do this first
        % start off with a weight of 0 at each vertex
        w = zeros(size(msh.vertex.origin,1),1);
        % Then loop over mesh vertices
        for ii = 1:length(msh_indices1)
            w(msh_indices1(ii)) = w(msh_indices1(ii))+weights(ii);
            w(msh_indices2(ii)) = w(msh_indices2(ii))+weights(ii);
        end
        
        % If the color range is not defined then go 1 to 90% of max
        if notDefined('crange')
            crange = [1 prctile(w(w>=1),90)];
        elseif length(crange) == 1
            crange = [crange(1) prctile(w(w>=1),90)];
        end
        
        % Dilate the coloring to adjacent vertices if desired
        if dilate == 1
            % Find faces that touch one of the colored vertices
            msh_faces = sum(ismember(msh.face.origin, [msh_indices1 msh_indices2]),2)>0;
            msh_indicesNew = msh.face.origin(msh_faces,:);
            msh_indicesNew = msh_indicesNew(:);
            % Confine this to only new mesh vertices that do not already have a fiber
            % endpoint
            msh_indicesNew = unique(msh_indicesNew(~ismember(msh_indicesNew,[msh_indices1 msh_indices2])));
            % Give these vertices a value as if one fiber is in them if there was zero
            % fibers. This effectively gives some smoothing
            for ii = 1:length(msh_indicesNew)
                % Find the adjacent vertices
                tmp = msh.face.origin(sum(msh.face.origin == msh_indicesNew(ii),2)>0,:);
                % And give our vertex in question the average, non-zero, weight of the
                % adjacent vertices
                wtmp = w(tmp(:));
                w(msh_indicesNew(ii)) = mean(wtmp(wtmp>0));
            end
        end
        
    case {'distance' 'dist' 'dweight' 'weighteddistance' 'weighteddist'}
        %% Find the distance from each endpoint to each mesh vertex
        
        % pull out the mesh vertices as a single to save space
        v_origin = single(msh.vertex.origin);
        % Convert start and endpoints to a single
        sp = single(sp); ep = single(ep);
        
        tic
        % First allocate space for a matrix of distances between each mesh
        % vertex (rows) and each fiber endpoint (columns). The 3rd
        % dimension is to consider start and endpointts separately. This
        % takes ~8 seconds for 3,000 fibers
        d2 = single(zeros([size(v_origin,1), size(sp,2), 2]));
        toc
        
        tic
        % Compute distance squared for each endpoint vertex combination. This takes
        % ~40 seconds.
        for ii = 1:length(sp)
            % for startpoints
            d2(:,ii,1) = 1./sum(bsxfun(@minus, v_origin, sp(:,ii)').^2,2);
            % for endpoints
            d2(:,ii,2) = 1./sum(bsxfun(@minus, v_origin, ep(:,ii)').^2,2);
        end
        toc
        
        tic
        % Threshold at 4mm. This takes rediculously long!
        d2(d2<1/16) = 0;
        toc
        
        % Normalize each column to sum to the prespecified weight for that
        % fiber. These weights will default to 1
        
end

%% Get the rgb color for each mesh vertex
rgb = vals2colormap(w, colors, crange);

% Color current mesh vertices for all vertices where the weight is above
% the defined threshold
msh.tr.FaceVertexCData(w>=crange(1),:) = rgb(w>=crange(1),:);

return



%% Notes of how to handle other vertex mappings

% First we need to check to see if the vertices are the same as the
% original ones
if  strcmp(msh.map2origin.(msh.vertex.current),'origin')
    % Combine colors
    facecolors = alpha.*repmat(color,length(msh_indices),1) + (1-alpha).*msh.tr.FaceVertexCData(msh_indices,:);
    msh.tr.FaceVertexCData(msh_indices,:) = facecolors;
else
    % Find the vertex to origin mapping
    new_indices = find(ismember(msh.map2origin.(msh.vertex.current), msh_indices));
    facecolors = alpha.*repmat(color,length(new_indices),1) + (1-alpha).*msh.tr.FaceVertexCData(new_indices,:);
    msh.tr.FaceVertexCData(new_indices,:) = facecolors;
end