function val = AFQ_meshGet(msh, param, varargin)
% Get parameters from an AFQ msh, mesh structure

param = mrvParamFormat(param);

switch(param)
    case {'triangles' 'tr'}
        val = msh.tr;
    case 'facevertexcdata'
        val = msh.tr.FaceVertexCData;
    case {'vertexorigin' 'originalvertices'}
        val = msh.vertex.origin;
    case 'basecolor'
        val = msh.colors.base;
    otherwise
        fprintf('\nNot an AFQ msh, mesh parameter');
       
end