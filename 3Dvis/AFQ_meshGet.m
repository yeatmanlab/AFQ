function val = AFQ_meshGet(msh, param, varargin)
% Get parameters from an AFQ msh, mesh structure


switch(param)
    case {'triangles' 'tr'}
        val = msh.tr;
    case 'FaceVertexCData'
        val = msh.tr.FaceVertexCData;
    case {'vertexorigin' 'originalvertices'}
        val = msh.vertex.origin;
    case 'basecolor'
        val = msh.colors.base;
    otherwise
        fprint('\nNot an AFQ msh, mesh parameter')
       
end