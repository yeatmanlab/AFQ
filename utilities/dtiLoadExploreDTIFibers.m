function fg = dtiLoadExploreDTIFibers(fgPath,dt)
% Load a fiber group from ExploreDTI
%
% fg = dtiLoadExploreDTIFibers(fgPath,dt)
%
% Inputs:
%
% fgPath = Path to ExploreDTI fiber group
% dt     = dt file converted from ExplorDTI. See dtiMakeDt6FromExploreDTI

load(fgPath)

for ii = 1:length(Tracts)
   fg.fibers{ii,1} = Tracts{ii}';
end
fg = dtiXformFiberCoords(fg,dt.xformToAcpc)