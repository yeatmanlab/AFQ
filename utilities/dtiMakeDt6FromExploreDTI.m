function dt = dtiMakeDt6FromExploreDTI(matfile, savefile)
% Make a dt6 from an explore DTI file. THIS FUNCTION IS STILL BEING
% DEVELOPED AND HAS BUGS!!
%
% t = dtiMakeDt6FromExploreDTI(matfile, savefile)

if ~exist('savefile') || isempty(savefile)
    savefile = 0;
end
% Load explore dti file
load(matfile)

% Make dt6
dt6 = zeros(size(DT{1},1),size(DT{1},2),size(DT{1},3),6,1);
dt6(:,:,:,:,1) = cat(4,DT{:});
dt6(isnan(dt6))=0;

dt.datafile = matfile;
dt.dt6 = dt6;
dt.b0 = DWIB0;
dt.xformToAcpc = diag([VDims, 1]);
dt.xformToAcpc(1:3,4) = -MDims./2;
dt.bb = [-80 80 -120; 90 -60 90];

if savefile == 1
    save dt6 dt
end