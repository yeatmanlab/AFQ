function TractProfile = AFQ_CreateTractProfile(varargin)
% Create the AFQ Tract Profile structure
%
% TractProfile = AFQ_CreateTractProfile(varargin)

% Check arguments
n = strcmpi('name',varargin);
if sum(n) > 0
    name = varargin{find(n)+1};
else
    name = {};
end

%% Create the afq tract profile structure
TractProfile.name = 'name';
TractProfile.vals   = struct('fa',[],'md',[],'rd',[],'ad',[],'cl',[]);
TractProfile.xform  = struct('acpc2mni',[],'mni2acpc',[],'acpc2dt6',[],'dt62acpc',[]);
TractProfile.coords = struct('acpc',[],'mni',[],'dt6',[]);


