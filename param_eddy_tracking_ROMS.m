% param_eddy_tracking_ROMS.m
%
%   param_eddy_tracking sets user defined paths and parameters for a
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
% Fixed parmeter:
%   - res: factor of interpolation of the fields (velocity and ssh)
%       to compute center detection. 2 or 3 seems reasonable in term of
%       computation time for 1/8° AVISO fields.
%   - K: LNAM(LOW<0) threshold to detect the potential eddy centers
%   - b: parameter for the computation of LNAM and Local Okubo-Weiss
%       (number of grid points in one of the 4 directions to define
%       the length of the box area used normalize Angular Momentum
%       and Okubo-Weiss fields)
%   - bx: number of grid points to define the initial area to scan
%       streamlines
%
%-------------------------
% IMPORTANT - Input file requirements:
%
% All the variables are read from netcdf file.
% The package requires 3 different input files:
% 1) nc_dim with variables x(j,i),y(j,i),mask(j,i)
% 2) nc_u with variable ssu(t,j,i), day(t)
% 3) nc_v with variable ssv(t,j,i), day(t)
%
% t is the temporal dimension (number of time steps)
% y is the zonal dimension (number of grid points along latitude)
% x is the meridional dimension (number of grid points along longitude)
%
% The grid is assumed to be rectangular, with orientation N-S and E-W. 
% Grid indices correspond to geography, so that point (1,1) represents the
%      south-western corner of the domain.
% Latitudinal and longitudinal grid spacing can vary within the grid domain.
%
%-------------------------
%   Ver Apr 2015 Briac Le Vu
%-------------------------
%
%=========================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User modification ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Experiment setings

% name for the experiment
name = 'Vel_Exp105_ROMS';

% set the paths
global path_in
global path_out
path_in = ['/home/blevu/DATA/ROMS/',name,'/'];
%path_out = ['/home/blevu/Resultats/ROMS/',name,'/tests/'];
path_out = ['/home/blevu/Resultats/ROMS/',name,'/tests_new/'];

% use to diferenciate source field of surface height (adt, ssh, psi,...)
global sshname
sshname = []; % ex: adt_

% use to submit parallel computation
global runname
runname = []; % ex: 1

% input data file absolute name
mat = [path_in,sshname,'velfield_day'];

% mesh grid  size (xd), deformation radius (Rd)
Dx = 2; % km
Rd = 5.9; % km

% rotation period (T) per day and time step in days (dps)
T = 3600*24; % day period in seconds
dps = 1; % 24h time step

% Calculate coriolis parameter
f = 4*pi/T; % in s-1

% degradation factor to test the algorithm
if ~exist('deg','var')
    deg = 1; % from 1 (default) to >10 in some experiment
end

%% Experiment option

% grid type
global grid_ll
grid_ll = 0;
        % 0 : spatial grid in cartesian coordinates (x,y)
        % 1 : spatial grid in earth coordinates (lon,lat)

% choose the field use as streamlines
global type_detection
type_detection = 1;
        % 1 : using velocity fields
		% 2 : using ssh
		% 3 : using both velocity fields and ssh, 
		%     and keep max velocity along the eddy contour

% if you want extended diags directly computed
global extended_diags
extended_diags = 1;
        % 0 : not computed
        % 1 : computed as the same time as eddy detection
		% 2 : computed afterward  

% save streamlines at days daystreamfunction and profil as well
global streamlines
streamlines = 0;
global daystreamfunction
daystreamfunction = 1:170;

% in case of periodic grid along x boundaries
global periodic
periodic = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of user modification ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Algortihm parametrisation from the test with AVISO, ROMS and PIV

% Resolution parameters:
%----------------------------------------------
% gama is resolution coefficient which is the number of pixels per Rd.
% After test gama>3 is required to get the max number of eddies.
gama = Rd / (Dx*deg); % around 1 for AVISO 1/8

% res is an integer and used to improve the precision of centers detection
% close to 3 pixels per Rd. res can goes up to 3
res = max(1,min(3,round(3/gama))); % [1 - 3]

% Detection parameters:
%----------------------------------------------
% K is LNAM(LOW<0) threshold to fixed contours where to detect potential center
% (one per contour). After test: no sensibility to K between 0.2 and 0.7
% but 0.7 give better time performance
K = 0.7; % [0 - 1]

% b is half length of the box in pixels used to normalise the LNAM and LOW.
% After test the optimal length of the box ( Lb = 2b*xd*deg )
% is fixed to one and half the size of Rd (Lb/Rd=1.2).
b = max(1,round((1.2*gama)/2));

% Rb (=Lb/Rd) is to check that the b setting and the gama are optimal.
% !!! Because optimal b is directly linked to gama, you start missing
% smaller eddies when gama is below 2 even at b=1
% (e.g. for AVISO 1/8 (gama~0.8) we have Rb~2.5) !!!
Rb = 2*b / gama;

% bx is half length of the box used in pixels to find contour streamlines
% in pixels bx = 10 times the Rd is enough to start testing any eddies
bx = max(1,round(10*gama));

% H is number of lines scanned in bx during detection of potential centers
% (centers0 to centers) and shapes (centers to centers2). After test 100
% seems a good compromise. This parameter is directly link to the time
% computation!
global H
H = 200; % number of line scanned

% Tracking parameters:
%----------------------------------------------
% parameter for the eddy tracking.
Dt = 7; % in days


