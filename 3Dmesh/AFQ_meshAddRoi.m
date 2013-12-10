function msh = AFQ_meshAddRoi(msh, roiPath, color, dilate, alpha)
% Color mesh vertices based on a region of interest

if notDefined('dilate')
    dilate = 0;
end
if notDefined('alpha')
    alpha = 1;
end
% Convert the nifti image to a set of coordinates
if ischar(roiPath)
    roi = dtiRoiFromNifti(roiPath, [],[],'mat',[],false);
    % Remove file extension and path to get the name of the image
    roiIm = readFileNifti(roiPath);
    [~,valname] = fileparts(roiIm.fname);
else
    roi = roiPath;
    valname = roi.name;
end
% Remove a secondary extension from the valname if there is one
valname = prefix(valname);
% remove any characters that are not allowed for field names
rmchar = {' ','_','1','2','3','4','5','6','7','8','9','0'};
for ch = 1:length(rmchar)
    if valname(1)==rmchar{ch}
        valname = horzcat('x',valname);
        continue
    end
end
valname(strfind(valname,' ')) = '_';
valname(strfind(valname,'-')) = '_';
% Find the closes mesh vertex to each coordinate
msh_indices = nearpoints(roi.coords', msh.vertex.origin');
% Dilate the roi to neighboring vertices
if dilate > 0
    for ii = 1:dilate
        % Find faces that touch one of the roi indices
        msh_faces = sum(ismember(msh.face.origin, msh_indices),2)>0;
        msh_indicesNew = msh.face.origin(msh_faces,:);
        % Add the vertices connected by this face to the roi
        msh_indices = unique(horzcat(msh_indices, msh_indicesNew(:)'));
    end
end
% Save these indices into a field titled based on the image name
msh.roi.(valname) = msh_indices;
msh.roi.show{end+1} = valname;

% Color current mesh vertices
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