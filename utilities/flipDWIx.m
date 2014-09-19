function [imPath,bvecsPath] =flipDWIx(imPath,bvecsPath)
% Flip the image and gradients over the x-axis and reset header
%
% This funciton is a hack to deal with Phillips data. For whatever reason
% the nifti's I got from the scanner have a flip in the header that causes
% problems somewhere in AFQ. This function deals with that

im = readFileNifti(imPath);
bvecs = dlmread(bvecsPath);

a = eye(4);
a(1,1) = -1;
for ii = 1:im.dim(4)
    im.data(:,:,:,ii) = affine(im.data(:,:,:,ii),a,[],0);
end

im.qto_xyz = a*im.qto_xyz;
im.qto_ijk = im.qto_ijk*a;
im.sto_ijk = im.sto_ijk*a;
im.sto_xyz = a*im.sto_xyz;
%im.qoffset_x = im.qoffset_x.*-1;
imPath = [prefix(prefix(imPath)) '_xflip.nii.gz'];
im.fname = imPath;
bvecsPath = [prefix(prefix(bvecsPath)) '_xflip.txt'];

writeFileNifti(im)

bvecs=mrAnatXformCoords(a,bvecs)';

dlmwrite(bvecsPath,bvecs)

