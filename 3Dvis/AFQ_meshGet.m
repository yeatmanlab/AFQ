function val = AFQ_meshGet(msh, param, varargin)
% Get parameters from an AFQ msh, mesh structure
%
% val = AFQ_meshGet(msh, param, varargin)
%
% AFQ_meshGet will get various parameters and values from the AFQ msh mesh
% structure. Examples:
%
% val = AFQ_meshGet('triangles') -- Get the the current faces, vertices and
% colors and return them in a structure that can be rendered with the
% MATLAB patch function.
%
% val = AFQ_meshGet('vertexcolors') -- Get the rgb values of the current
% mesh vertices.
%
% val = AFQ_meshGet('vertexorigin') -- Get the original vertices, before
% any smoothing. These correspond to acpc coordinates of the original
% segmentaion image.
%
% val = AFQ_meshGet('basecolor' 
% val = AFQ_meshGet('colors', colorname
% val = AFQ_meshGet('currentvertexname'
%
% Copyright Jason D. Yeatman, November 2012
param = mrvParamFormat(param);

switch(param)
    case {'triangles' 'tr'}
        val = msh.tr;
    case {'vertexcolors' 'facevertexcdata'}
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