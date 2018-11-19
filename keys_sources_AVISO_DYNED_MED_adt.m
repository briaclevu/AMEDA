% keys_sources_AVISO.m
%
%   keys_sources sets user defined paths and user keys
%
% Paths:
%   - path_in: directory containing the input files;
%              (default is '..\Data\')
%   - path_out: directory for the output files;
%              (default is '..\Tracks\')
%   - nc_u: full name of the netcdf file with the zonal component of
%           velocity (ssu) and the time index (step)
%   - nc_v: full name of the netcdf file with the meridional component of
%           velocity (ssv) and the time index (step)
%   - nc_dim: full name of the netcdf file with the domain coordinates 
%            (longitude and latitude) and the velocity mask (land-points=0;
%             ocean-points=1)
%
% User option keys:
%   - type_detection: flag to choose the field use as streamlines
%           1 : using velocity fields
%           2 : using ssh
%           3 : using both velocity fields and ssh for eddy contour, 
%               and keep the contour of the last one if it existS.
%   - extended_diags: flag to have extra diags concerning eddies directly
%       computed (like ellipse features or vorticity for each eddy)
%           0 : not computed
%           1 : computed as the same time as eddy detection
%           2 : computed afterward
%   - streamlines and daystreamfunction: save streamlines at steps 
%       'daystreamfunction' and profils of streamlines scanned as well
%       (1:stepF by default )
%   - periodic: flag to activate options for East-West periodic
%               (e.g. global fields or idealized simulations) domains.
%               IMPORTANT: first and last columns of the domain must be
%                          identical for this to work properly!!!!!
%
%-------------------------
% IMPORTANT - Input file requirements:
%
% All the variables are read from netcdf file.
% The package requires 3 different input files:
% 1) nc_dim with variables x(j,i),y(j,i) and mask(y,x)
% 2) nc_u with variable ssu(t,j,i) in m/s, step(t)
% 3) nc_v with variable ssv(t,j,i) in m/s, step(t)
% you can have also the sea level in an other file
% 4) nc_ssh with variable ssh(t,j,i) in m, step(t)
%
% t is the temporal dimension (number of time steps)
% j is the zonal dimension (number of grid points along y or latitude)
% i is the meridional dimension (number of grid points along x or longitude)
%
% The grid is assumed to be rectangular, with orientation N-S and E-W. 
% Grid indices correspond to geography, so that point (1,1) represents the
%      south-western corner of the domain.
% Latitudinal and longitudinal grid spacing can vary within the grid domain.
%
%-------------------------
%   Ver Jun 2018 Briac and Romain Pennel
%   Ver Apr 2015 Briac Le Vu
%   Ver. 2.1 Oct.2012
%   Ver. 2.0 Jan.2012
%   Ver. 1.3 Apr.2011
%   Ver. 1.2 May.2010
%   Ver. 1.1 Dec.2009
%   Authors: Francesco Nencioli, francesco.nencioli@univ-amu.fr
%            Charles Dong, cdong@atmos.ucla.edu
%-------------------------
%
%=========================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User modification ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Experiment setings
%----------------------------------------------

% source of the input
source='AVISO';

% postfix name of the data
postname = '_2000_2017';

% set name of the domain
config='DYNED_MED';

% use to diferenciate source field of surface height (adt, ssh, psi,...)
sshtype='adt'; % adt or sla

% use to differenciate job
runname = '';

% set the paths
path_in = ['/home/blevu/DATA/',source,'/',config,'/'];
path_result = ['/home/blevu/Resultats/',source,'/',config,'/'];
path_out = [path_result,sshtype,postname,'/',runname,'/'];
path_tracks1 = ['/home/blevu/DATA/',source,'/dt/adt/tracks/'];
path_tracks2 = ['/home/blevu/DATA/',source,'/nrt/adt/tracks/'];
path_data = '/home/blevu/DATA/CORIOLIS/';
path_rossby = '/home/blevu/MATLAB/Rossby_radius/';

disp(['Compute from ',path_in])
disp([' to ',path_out])

if exist(path_out,'file')==0
    system(['mkdir ',path_result]);
    system(['mkdir ',path_result,sshtype,postname,'/']);
    system(['mkdir ',path_result,sshtype,postname,'/',runname,'/']);
    system(['mkdir ',path_out]);
end

% add path_out path
addpath(path_out)

% input data file absolute name
nc_dim = [path_in,'lon_lat_',sshtype,'_',config,'.nc'];
nc_u = [path_in,'ssu_',sshtype,'_',config,postname,'.nc'];
nc_v = [path_in,'ssv_',sshtype,'_',config,postname,'.nc'];
nc_ssh = [path_in,'ssh_',sshtype,'_',config,postname,'.nc'];

% Rossby deformation radius file
mat_Rd = [path_rossby,'Rossby_radius']; % for Med
name_Rd = 'Rd_baroc1_extra';

% searched eddies typical radius
% eddies size will be restricted to 5 times this radius
% eddies smaller than 1/4 this radius will be smoothed
Rd_typ = 12; % km

% minimal size for rmax to be reasonably detected
nRmin = 1/2; % half of the native Dx grid size

% variable names (could be automatised)
y_name = 'lat';
x_name = 'lon';
m_name = 'mask';
u_name = 'u';
v_name = 'v';
s_name = 'ssh';

% duration experiment (should read a 'time' variable)
if ~exist('stepF','var')
    u0 = squeeze(ncread(nc_u,u_name,[1 1 1],[1 1 Inf]));
    stepF = length(u0);
    clear u0
end

disp([' ',num2str(stepF),' time steps'])
disp(' ')

% rotation period (T) per day and time step in days (dps)
T = 3600*24; % day period in seconds
dps = 1; % 24h time step
level = 1; % vertical level (used only if any)

% degradation factor to test the algorithm
if ~exist('deg','var')
    deg = 1; % from 1 (default) to >10 in some experiment
end

%% Experiment option keys
%----------------------------------------------

% grid type
grid_ll = 1;
        % 0 : spatial grid in cartesian coordinates (x,y)
        % 1 : spatial grid in earth coordinates (lon,lat)

% grid regular or not (like arakawa in NEMO)
grid_reg = 1;
        % 0 irregular 
        % 1 regular 

% choose the field use as streamlines
type_detection = 1;
        % 1 : using velocity fields
        % 2 : using ssh
        % 3 : using both velocity fields and ssh, 
        %     and keep max velocity along the eddy contour

% if you want extended diags directly computed
extended_diags = 1;
        % 0 : not computed
        % 1 : computed as the same time as eddy detection
        % 2 : computed afterward  

% save streamlines at days daystreamfunction and profil as well
streamlines = 1;
daystreamfunction = 1:stepF;

% in case of periodic grid along x boundaries
periodic = 0;

% to keep firts and last detection after the tracking in NRT configuration
nrt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of user modification ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
