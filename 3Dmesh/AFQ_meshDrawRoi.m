function [coords, indices, bin, msh] = AFQ_meshDrawRoi(fh, msh, dilate)
%
%
% msh = AFQ_meshCreate;
% [msh.tr.vertices, msh.tr.faces] = read_surf('/mnt/diskArray/projects/freesurfer/fsaverage/surf/lh.white');
% msh.tr.faces = msh.tr.faces + 1;
% msh.tr.FaceVertexCData = repmat([.8 .7 .6],size(msh.tr.vertices,1),1);
% msh.tr = smoothpatch(msh.tr,[],20);
% [coords, indices, bin, msh] = AFQ_meshDrawRoi([],msh, 7)
% AFQ_RenderCorticalSurface(msh);
if ~exist('fh', 'var') || isempty(fh)
    [~,~,l]=AFQ_RenderCorticalSurface(msh);
    view(180,-90); camlight(l,'infinite');
    fh = gcf;
    fh = fh.Number;
end
datacursormode on
dcmObj = datacursormode(fh);
set(dcmObj,'SnapToDataVertex','on','Enable','on');
keypress = 0; ii = 0;
hold('on');
while keypress==0
    ii = ii+1;
    [keypress] = waitforbuttonpress;
    point = getCursorInfo(dcmObj);
    coords(ii,:) = point.Position
    plot3(coords(ii,1),coords(ii,2),coords(ii,3),'ko','markerfacecolor','k');
end
hold('off');

% get index of the coords
[~,indices] =   ismember(coords, msh.tr.vertices, 'rows');

if exist('dilate' ,'var') && ~isempty(dilate) && dilate>0
    for d = 1:dilate
        for ii = 1:length(indices)
            r = msh.tr.faces==indices(ii);
            tmp = msh.tr.faces(find(sum(r,2)>0),:);
            indices = unique(vertcat(indices,tmp(:)));
        end 
    end
end

% close the figure
close(fh)
msh.tr.FaceVertexCData(indices,:) = repmat([1 0 0],length(indices),1);

% Make a logical vector out of the indices
bin = zeros(size(msh.tr.vertices,1),1);
bin(indices) = 1;
