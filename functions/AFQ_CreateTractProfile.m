function TractProfile = AFQ_CreateTractProfile(varargin)
% Create the AFQ Tract Profile structure
%
% TractProfile = AFQ_CreateTractProfile(varargin)
    
%% Create the afq tract profile structure
TractProfile.name    = 'name';
TractProfile.nfibers = [];
TractProfile.vals    = struct('fa',[],'md',[],'rd',[],'ad',[],'cl',[],'planarity',[],'sphericity',[]);
TractProfile.xform   = struct('acpc2mni',[],'mni2acpc',[],'acpc2dt6',[],'dt62acpc',[]);
TractProfile.coords  = struct('acpc',[],'mni',[],'dt6',[]);
TractProfile.fibercov     = [];
TractProfile.fiberCurvature = [];
TractProfile.fiberTorsion = [];
TractProfile.fiberVolume = [];
TractProfile.fiberCovVolume = [];
%% Set fields
TractProfile = afqVarargin(TractProfile,varargin);
%% Check if the old super fiber structure was passed in
if nargin > 0 && sum(strcmpi('superfiber',varargin)) == 1
    idx = find(strcmpi('superfiber',varargin));
    TractProfile.fibercov = varargin{idx+1}.fibervarcovs{1};
    TractProfile.coords.acpc = varargin{idx+1}.fibers{1};
    TractProfile.nfibers = varargin{idx+1}.n;
end
    

