function msh = AFQ_meshSet(msh, param, varargin)
% Set fields in the AFQ msh mesh structure
%
% msh = AFQ_meshSet(msh, 'filtered', 'box5', tr)
% msh = AFQ_meshSet(msh, 'vertices', 'filtered box 5')
% msh = AFQ_meshSet(msh, 'vertex','name', vertices)
% msh = AFQ_meshSet(msh, 'vertex','name', vertices, faces)
% msh = AFQ_meshSet(msh, 'vertex','name', vertices, faces, map2origin)
switch(param)
    case('vertices')
        % Set which vertices to render
        if nargin < 3, fprintf('Please supply vertex name'); end
        % If that vertex type is a field then set it to tr.vertices. If
        % it's not a field but can be computed then compute it and add it to
        % tr.vertices
        if isfield(msh.vertex, varargin{1}) && isempty(strfind(varargin{1}, 'filtered')) 
            msh.tr.vertices = msh.vertex.(varargin{1});
            % if the faces are not set properly then set them
            if isnumeric(msh.face.(varargin{1}))
                msh.tr.faces = msh.face.(varargin{1});
            elseif ischar(msh.face.(varargin{1})) && ...
                    strcmp(msh.face.(varargin{1}), 'origin')
                msh.tr.faces = msh.face.origin;
            end
            
        elseif  isfield(msh.vertex, varargin{1}) && strfind(varargin{1}, 'filtered') == 1
            % The filtered vertices are stored slightly differntly
            msh.tr.vertices = msh.vertex.(varargin{1}).vertices;
            msh.tr.faces     = msh.vertex.(varargin{1}).faces;

        elseif strfind(varargin{1},'smooth') == 1
            % If the desired vertices do not exist then compute them
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
    case {'filter' 'filtered'}
        fname = ['filtered' varargin{1}];
        msh.vertex.(fname) = varargin{2};
    case 'vertex'
        if nargin < 4
            fprintf('Correct call: msh = AFQ_meshSet(msh,''vertex'',''name'',vertices, faces, map2origin)\n')
        elseif nargin == 4
            msh.vertex.(varargin{1}) = varargin{2};
            % if no faces were passed in then assume that the faces are the
            % same the original faces
            msh.face.(varargin{1})       = 'origin';
            msh.map2origin.(varargin{1}) = 'origin';
        elseif nargin == 5
            msh.vertex.(varargin{1}) = varargin{2};
            msh.face.(varargin{1})   = varargin{3};
            if ~strcmp(varargin{1},'origin')
                fprintf('no data for ''map2origin''. Mesh surface may not be colored properly\n');
            end
        elseif nargin == 6
            msh.vertex.(varargin{1})     = varargin{2};
            msh.face.(varargin{1})       = varargin{3};
            msh.map2origin.(varargin{1}) = varargin{4};
        end
end