function msh = AFQ_meshCut(msh, plane)
%
%
% Inputs:
% msh
% plane   - 2x3 matrix where each column is x,y,z and the row is the clip
%           range

cp = diff(plane);

% X clip
if cp(1) >0
    % These are the vertices to remove. All the ones that are greater than
    % the first plane but less than the second plane
    vInds = msh.tr.vertices(:,1) > plane(1,1) & msh.tr.vertices(:,1) < plane(2,1);
    % Remove vertices and corresponding colors within the clip range
    msh.tr.vertices = msh.tr.vertices(~vInds,:);
    msh.tr.FaceVertexCData = msh.tr.FaceVertexCData(~vInds,:);
    % Find faces that include a removed vertex
    fInds = [];
    fInds(:,1) = ismember(msh.tr.faces(:,1),find(vInds));
    fInds(:,2) = ismember(msh.tr.faces(:,2),find(vInds));
    fInds(:,3) = ismember(msh.tr.faces(:,3),find(vInds));
    % Turn into a logical vector
    Fkeep = sum(fInds,2) == 0;
    % Remove any face that goes to a removed vertex
    msh.tr.faces = msh.tr.faces(Fkeep,:);
    % Now this part is a bit tricky. We need to recompute which rows of
    % vertices each face corresponds to since we have removed a bunch of
    % rows
    offset = cumsum(vInds);
    for ii = 1:numel(msh.tr.faces)
        msh.tr.faces(ii) = msh.tr.faces(ii) - offset(msh.tr.faces(ii));
    end
end
% Y clip
if cp(2) >0
    % These are the vertices to remove. All the ones that are greater than
    % the first plane but less than the second plane
    vInds = msh.tr.vertices(:,2) > plane(1,2) & msh.tr.vertices(:,2) < plane(2,2);
    % Remove vertices and corresponding colors within the clip range
    msh.tr.vertices = msh.tr.vertices(~vInds,:);
    msh.tr.FaceVertexCData = msh.tr.FaceVertexCData(~vInds,:);
    % Find faces that include a removed vertex
    fInds = [];
    fInds(:,1) = ismember(msh.tr.faces(:,1),find(vInds));
    fInds(:,2) = ismember(msh.tr.faces(:,2),find(vInds));
    fInds(:,3) = ismember(msh.tr.faces(:,3),find(vInds));
    % Turn into a logical vector
    Fkeep = sum(fInds,2) == 0;
    % Remove any face that goes to a removed vertex
    msh.tr.faces = msh.tr.faces(Fkeep,:);
    % Now this part is a bit tricky. We need to recompute which rows of
    % vertices each face corresponds to since we have removed a bunch of
    % rows
    offset = cumsum(vInds);
    for ii = 1:numel(msh.tr.faces)
        msh.tr.faces(ii) = msh.tr.faces(ii) - offset(msh.tr.faces(ii));
    end
end