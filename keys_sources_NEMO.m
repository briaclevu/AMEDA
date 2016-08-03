% keys_sources_NEMO.m
%
%   keys_sources sets user defined paths and parameters for a
%   degradation coefficient of 'deg' which goes from 1 (default)
%   to >10 in some experiment.
%
% Paths:
%   - path_in: directory containing the input files;
%              (default is '..\Data\')
%   - path_out: directory for the output files;
%              (default is '..\Tracks\')
%   - nc_u: full name of the netcdf file with the zonal component of
%           velocity (ssu) and the time index (day)
%   - nc_v: full name of the netcdf file with the meridional component of
%           velocity (ssv) and the time index (day)
%   - nc_dim: full name of the netcdf file with the domain coordinates 
%            (longitude and latitude) and the velocity mask (land-points=0;
%             ocean-points=1)
%
% User parameters definition:
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
% 1) nc_dim with variables lon(j,i),lat(j,i) and mask(y,x)
% 2) nc_u with variable ssu(t,y,x) in m/s, day(t)
% 3) nc_v with variable ssv(t,y,x) in m/s, day(t)
% 4) nc_ssh with variable ssh(t,y,x) in m, day(t)
%
% t is the temporal dimension (number of time steps)
% j is the zonal dimension (number of grid points along latitude)
% i is the meridional dimension (number of grid points along longitude)
%
% The grid is assumed to be rectangular, with orientation N-S and E-W. 
% Grid indices correspond to geography, so that point (1,1) represents the
%      south-western corner of the domain.
% Latitudinal and longitudinal grid spacing can vary within the grid domain.
%-------------------------
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

% name for the experiment
name = 'FOR_12';

% use to diferenciate source field of surface height (adt, ssh, psi,...)
sshname = 'psi_'; % ex: adt_

% set the paths
path_ameda = ['/home/blevu/MATLAB/AMEDA_v2/'];
path_in = ['/home/blevu/DATA_tmp/NEMO/',name,'/'];
path_out = ['/home/blevu/Resultats/NEMO/',name,'/'];
%path_out = ['/data/BIGWORKSPACE/tarsouze/WMED36/detection/30m_v2/'];

% use to submit parallel computation
runname = 'jan_mar13'; % ex: 1,b2_K5,day10,...

% input data file absolute name
%nc_u=[path_in,'EXP2FL_3h_20130101_20130315_grid_U_zoom_30m_day.nc'];
%nc_v=[path_in,'EXP2FL_3h_20130101_20130315_grid_V_zoom_30m_day.nc'];
%nc_dim=[path_in,'mesh_mask_WMED36_zoom_30m.nc'];
%nc_ssh=[path_in,'EXP2FL_3h_20130101_20130315_ssh_zoom_day.nc'];
nc_u = [path_in,'EXP2FL_3h_20130101_20130315_grid_U_zoom_30m.nc'];
nc_v = [path_in,'EXP2FL_3h_20130101_20130315_grid_V_zoom_30m.nc'];
nc_dim = [path_in,'mesh_mask_WMED36_zoom_30m.nc'];
nc_ssh = [path_in,'EXP2FL_3h_20130101_20130315_ssh_zoom.nc'];

% variable names
y_name  = 'nav_lat';
x_name  = 'nav_lon';
m_name = 'tmask';
u_name = 'vozocrtx';
v_name = 'vomecrty';
s_name= 'sossheig';

% duration experiment
u0 = squeeze(ncread(nc_u,u_name,[1 1 1 1],[1 1 1 Inf]));
stepF = length(u0);
clear u0

% rotation period (T) and day per time step (dps)
T = 3600*24; % 24h
%dps = 1; % 1 is daily time step
dps = 1/8; % 1/8 is 3h time step

% degradation factor to test the algorithm
if ~exist('deg','var')
    deg = 1; % from 1 (default) to >10 in some experiment
end

%% Experiment option
%----------------------------------------------

% grid type
grid_ll = 1;
        % 0 : spatial grid in cartesian coordinates (x,y)
        % 1 : spatial grid in earth coordinates (lon,lat)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of user modification ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

