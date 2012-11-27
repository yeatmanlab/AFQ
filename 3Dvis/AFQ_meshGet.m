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
    case {'basecolor' 'base'}
        val = msh.colors.base;
    case {'colors' 'color'}
        if nargin == 3 && isfield(msh.colors,varargin{1}) 
            % Get the color
            val = msh.colors.(varargin{1});
            
            % If the color is a single rgb value the repeat it for each
            % vertex
            if size(val,1) == 1
                val = repmat(val,size(msh.vertex.(AFQ_meshGet(msh,'currentvertexname')),1),1);
            
            % If it is a matrix of rgb values but not the same size as the vertices then map the
            % correct colors with the map2origin field
            elseif size(val,1) ~= size(msh.tr.vertices,1)
                val = val(msh.map2origin.(AFQ_meshGet(msh,'currentvertexname')),:);
                if size(val) ~= size(msh.tr.vertices)
                    error('mismatch between size of vertices and coloring');
                end
            end
        else
            error('\nPlease provide the name of a color in msh.colors')
        end
        
    case 'currentvertexname'
        val = msh.vertex.current;
    otherwise
        error('\nNot an AFQ msh, mesh parameter');
       
end