function imOut = niftiChangeResolution(im, mmPerVox)

% First we build a bounding box in pixels (image space) containing the full
% volume of the second input image.
bbimg = [1 1 1; im.pixdim(1:3)./mmPerVox.*im.dim(1:3)];

% Second, we convert the bounding box from image space to mm space.
bbmm = mrAnatXformCoords(finalResNifti.qto_xyz,bbimg);
[newImg,xform,deformField] = mrAnatResliceSpm(im.data, im.qto_xyz, bb, mmPerVox, bSplineParams, showProgress)

% First we create a blank nifti image with the desired output resolution.
finalResNifti = niftiCreate;

finalResNifti.pixdim = mmvox;
finalResNifti.dim = round(im.pixdim./mmvox.*im.dim)
finalResNifti.qto_xyz 


resample_params = [7 7 7 0 0 0];
% Resample the image to match the new blank one
imOut = mrAnatResampleToNifti(im, finalResNifti,[],resample_params);
