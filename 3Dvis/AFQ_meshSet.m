function msh = AFQ_meshSet(msh, param, varargin)
% Set fields in the AFQ msh mesh structure
%
% msh = AFQ_meshSet(msh, 'filter', 'box5', tr)

switch(param)
    case('vertices')
        % Set which vertices to render
        if nargin < 3, fprintf('Please supply vertex name'); end
        % If that vertex type is a field then set it to tr.vertices. If
        % it's not a field but can be computed then compute it and add it to
        % tr.vertices
        if isfield(msh.vertex, varargin{1})
            msh.tr.vertices = msh.vertex.(varargin{1});
        elseif strfind(varargin{1},'smooth') == 1
            sIter = str2num(varargin{1}(7:end));
            vName = sprintf('smooth%d',sIter);
            % Set the vertices to the original ones so that they can be
            % smoothed
            msh = AFQ_meshSet(msh,'vertices','origin');
            % Smooth them
            tmp = smoothpatch(AFQ_meshGet(msh,'tr'), [], sIter);
            % Add them to the corresponding vertex field
            msh.vertex.(vName) = tmp.vertices;
            % Set these new vertices to the ones that will be rendered
            msh = AFQ_meshSet(msh,'vertices',vName);
        end
    case 'FaceVertexCData'
        % If a single color value is in the color field and no color name
        % was specified then replicate it for al the vertices
        if nargin < 3 && size(AFQ_meshGet(msh,'FaceVertexCData'),1) == 1
            msh.tr.FaceVertexCData = repmat(msh.tr.FaceVertexCData,size(msh.tr.vertices,1),1);
        else
            msh.tr.FaceVertexCData = msh.colors.(varargin{1});
        end
        % If there was only a single color provided then replicate it for
        % all vertices
        if size(AFQ_meshGet(msh,'FaceVertexCData'),1) == 1
            msh.tr.FaceVertexCData = repmat(AFQ_meshGet(msh,'FaceVertexCData'),size(msh.tr.vertices,1),1);
        end
    case {'vals' 'val'}
        if nargin < 4
            fprintf('Please provide the value name and associated vector of values\n')
        elseif ischar(varargin{1}) && isvector(varargin{2})
            msh.vals.(varargin{1}) = varargin{2};
        else
            fprintf('Correct call: msh = AFQ_meshSet(msh,''vals'',''valname'',values)\n');
        end
    case 'colors'
        if nargin < 4
            fprintf('Please provide the color name and associated vector of color values)\n')
        elseif ischar(varargin{1}) && ismatrix(varargin{2}) && size(varargin{2},2) == 3;
            msh.colors.(varargin{1}) = varargin{2};
        else
            fprintf('Correct call: msh = AFQ_meshSet(msh,''colors'',''valname'',values)\n');
        end
    case {'basecolor' 'color'}
        if length(varargin) == 1 && size(varargin{1},2) == 3
            msh.colors.base = varargin{1};
        else 
            fprintf('Correct call: msh = AFQ_meshSet(msh,''basecolor'',rgbValue)\n');
        end
    case 'filter'
        msh.filter.(varargin{1}) = varargin{2};
end