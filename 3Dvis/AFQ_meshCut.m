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
    vInds = msh.tr.vertices(:,1) > plane(1,1) & msh.tr.vertices(:,1) < plane(2,1);
    % Remove vertices and corresponding colors within the clip range
    msh.tr.vertices = msh.tr.vertices(~vInds,:);
    msh.tr.FaceVertexCData = msh.tr.FaceVertexCData(~vInds,:);
    % Find face corresponding to these vertices
    fInds(:,1) = ismember(msh.tr.faces(:,1),find(vInds));
    fInds(:,2) = ismember(msh.tr.faces(:,2),find(vInds));
    fInds(:,3) = ismember(msh.tr.faces(:,3),find(vInds));
    % Remove any face that goes to a removed vertex
    msh.tr.faces = msh.tr.faces(sum(fInds,2)==0);
end