function msh = AFQ_meshAddRoi(msh, roiPath, color)
% Color mesh vertices based on a region of interest

% Convert the nifti image to a set of coordinates
roi = dtiRoiFromNifti(roiPath, [],[],'mat',[],false);

% Find the closes mesh vertex to each coordinate
msh_indices = nearpoints(roi.coords', msh.vertex.origin');

% Remove file extension and path to get the name of the image
roiIm = readFileNifti(roiPath);
[~,valname] = fileparts(roiIm.fname);
% Remove a secondary extension if there is one
valname = prefix(valname);
% Save these indices into a field titled based on the image name
msh.roi.(valname) = msh_indices;
msh.roi.show{end+1} = valname;

% Color current mesh vertices
% First we need to check to see if the vertices are the same as the
% original ones
if  strcmp(msh.map2origin.(msh.vertex.current),'origin')
    msh.tr.FaceVertexCData(msh_indices,:) = repmat(color,length(msh_indices),1);
else
    % Find the vertex to origin mapping
    new_indices = find(ismember(msh.map2origin.(msh.vertex.current), msh_indices));
    msh.tr.FaceVertexCData(new_indices,:) = repmat(color,length(new_indices),1);
end