function mod_eddy_params(keys_sources,stepF)
%mod_eddy_params(keys_sources {,stepF})
%
%   mod_eddy_params loads user paths and keys defined in keys_sources.m
%   and sets 2D fields of parameters
%
% Computed parameters:
%   - b: parameter for the computation of LNAM and Local Okubo-Weiss
%       (number of grid points in one of the 4 directions to define
%       the length of the box area used normalize Angular Momentum
%       and Okubo-Weiss fields)
%   - bx: number of grid points to define the initial area to scan
%       streamlines
%   - Dx: Meshgrid size
%   - Rd: First Baroclinic Rossby Radius of Deformation
%   - gama: resolution coefficient which is the number of pixels per Rd.
%   - resol: factor of interpolation of the fields (velocity and ssh)
%       to compute center detection. 2 or 3 seems reasonable in term of
%       computation time for 1/8° AVISO fields.
%
% Fixed parameters
% (you can play with Dt and cut_off but don't touch too much others):
%   - K: LNAM(LOW<0) threshold to detect the potential eddy centers
%   - r: center amximale speed during tracking of an eddy
%   - Dt: delay in days for the tracking tolerance
%   - cut_off: in days for eddies duration filtration after the tracking
%           0 : use the turnover time from each eddies
%           1 : keep all the tracked eddies
%   - D_stp: in steps to define the averaged eddy features during the tracking
%   - N_scan: maximum number of candidat the assignment must resolve
%   - C_rad: number of the averaged radius for the search limit area
%   - H: number of scanned streamlines in the initial area
%       for each eddies (see 'box' parameter below)
%   - n_min: minimal number of point to defined a contour as streamlines
%   - vel_epsilon: minimal difference between 2 tests of the velocity
%       to admit an increase
%   - k_vel_decay: coefficient of velmax to detect a decrease
%   - R_lim: limite of the eddy size
%   - nrho_lim: limite of the negative curvature for a single eddy
%   - ds_max: maximal distance between 2 eddies centers to have a potential
%       interaction (merging or splitting)
%   - dc_max: maximal distance between 2 eddies contours
%
%----------------------------------------------
% This parametrisation come from the test with AVISO, ROMS and PIV
% see LE VU et al. paper
%
%-------------------------
%   June 2016 Briac Le Vu
%-------------------------
%
%=========================

%% Load keys and paths file
%----------------------------------------------
run(keys_sources)

%% Calculate the meshgrid size at (x,y) coordinates
%----------------------------------------------

% Read grid
lon0 = double(ncread(nc_dim,x_name))';
lat0 = double(ncread(nc_dim,y_name))';

% produce degraded field 
x = lon0(1:deg:end,1:deg:end);
y = lat0(1:deg:end,1:deg:end);

[N,M] = size(x);

% Meshgrid size
if grid_ll
    Dx = get_Dx_from_ll(x,y);
else
    Dx = ( abs( diff(x,1,2) ) + abs( diff(y,1,1) ) )/2;
end

% Calculate coriolis parameter
f = 4*pi/T; % in s-1

g=9.8;

%% fixed parameters
%----------------------------------------------

% Scanning parameters:
%----------------------------------------------
% H is number of contour scanned in box during detection of potential centers
% (centers0 to centers) and shapes (centers to centers2).
% After test [100-200] seems a good compromise.
% This parameter is directly link to the time computation!
H = 200; % number of contour scanned

% these contours must contain a minimal number of dots ([4-6])
n_min = 4; % minimal number of point to defined a contour as streamlines

% minimal number of point to calculate curvature along a segment at point i
Np = 3; % (segment i-Np:i+Np)

% minimal difference ([1-2]*1e-4) between 2 tests of the velocity
% to admit an increase
vel_epsil = 1e-4; % in m/s

% coefficient of velmax to detect a decrease ([0.95-0.99])
k_vel_decay = 0.95; % in term of velmax

% size limite in term of deformation radius ([4-7]Rd) for a single eddy
nR_lim = 5; % [4-7]Rd

% limite of the length of the contour with negative curvature for a single
% eddy from empirical determination based on ROMS, PIV and AVISO tests
nrho_lim = 0.2; % [1/5-1/2]*2*pi*Rmax of the equivalent circle

% Double eddy parameters
%----------------------------------------------
% maximal distance distances between 2 eddies centers to have a potential
% interaction (merging or splitting) are based on von Hardenberg et al. (2000).
% study and synthesis on vortex merging in barotropic flows
% (ie d=3.3Rmax the critical merging distance).

% maximal distance between 2 eddies centers in term of Rmax
dc_max = 3.5; 

% Tracking parameters:
%----------------------------------------------
% maximal eddy tracking distance between 2 time steps.
%!!! can be change to vmax for each eddy !!!
r = 6.5; % 6.5km/day

% maximal delay after its last detection for tracking an eddy [1-15]
% represent the temporal correlation depending on the coverage of an area
% 2 steps for model, 3-5 steps for imagery experiment.
%!!! in case of AVISO this will be ajusted with the error map of aviso !!!
Dt = 10; % in days

% minimal duration of recorded eddies [0-100]
% 0 for keeping only eddies longer than their turnover time
cut_off = 0; % in days

% number of steps to define the averaged eddy features during the tracking
D_stp = 4; % in steps

% number of candidats during the tracking (better to be high)
N_can = 20; % in steps

%% Calculate parameters needed for AMEDA at Dx and a given Rd
%----------------------------------------------

% Interpolate First Baroclinic Rossby Radius of Deformation
% (taken from 2D file Rossby_radius computed using Chelton et al. 1998)
%----------------------------------------------
load(mat_Rd)
Rd = interp2(lon_Rd,lat_Rd,Rd_baroc1_extra,x,y); % 10km in average AVISO 1/8

% Resolution parameters:
%----------------------------------------------
% gama is resolution coefficient which is the number of pixels per Rd.
% After test gama>3 is required to get the max number of eddies.
gama = Rd ./ Dx; % [0.1-1.5] for AVISO 1/8 (0.8 in average)

% resol is an integer used to improve the precision of centers detection
% close to 3 pixels per Rd. resol can goes up to 3
resol = max(1,min(3,round(3/min(gama(:))))); % [1 - 3]

% Detection parameters (cf. Le Vu et al. 2016 for test results):
%----------------------------------------------
% K is LNAM(LOW<0) threshold to fixed contours where to detect potential center
% (one per contour). After test: no sensibility to K between 0.2 and 0.7
% but 0.7 give better time performance
K = 0.7; % [.4-.8]

% b is half length of the box in pixels used to normalise the LNAM and LOW.
% After test the optimal length of the box ( Lb = 2b*Dx*deg )
% is fixed to 1.2 the size of Rd (Lb/Rd=1.2).
b = max(1,round((1.2*gama)/2)); % always 1 for AVISO 1/8

% Rb (=Lb/Rd) is to check that the b setting and the gama are optimal.
% !!! Because optimal b is directly linked to gama, you start missing
% smaller eddies when gama is below 2 even at b=1
% (e.g. for AVISO 1/8 (gama~0.8) we have Rb~2.5) !!!
Rb = 2*b ./ gama; % [1-100] for AVISO 1/8 (2.5 in average)

% box is half length of the box used in pixels to find contour streamlines
% in pixels box = 10 times the Rd is enough to start testing any eddies
bx = max(1,round(2*nR_lim*gama)); % [1-14] for AVISO 1/8 (8 in average)

%% interpolated parameters
%----------------------------------------------

if resol==1

    bi = b;
    bxi = bx;
    Dxi = Dx;
    Rdi = Rd;
    gamai = gama;

else

    % size for the interpolated grid
    Ni = resol*(N-1)+1; % new size in lat
    Mi = resol*(M-1)+1; % new size in lon

    % elemental spacing for the interpolated grid with regular grid
    dy = diff(y(1:2,1))/resol;
    dx = diff(x(1,1:2))/resol;

    % interpolated grid
    [xi,yi] = meshgrid([0:Mi-1]*dx+min(x(:)),[0:Ni-1]*dy+min(y(:)));

    % Compute interpolated b and bx
    bi = round(interp2(x,y,b,xi,yi))*resol;
    bxi = round(interp2(x,y,bx,xi,yi))*resol;

    % Dx, Rd, gama on interpolated grid
    Dxi = get_Dx_from_ll(xi,yi);
    Rdi = interp2(lon_Rd,lat_Rd,Rd_baroc1_extra,xi,yi); % 10km in average AVISO 1/8
    gamai = Rdi ./ Dxi / resol;

end
    
%% Save parameters and paths
%----------------------------------------------
save([path_out,'param_eddy_tracking'],...
    'postname','domain','sshtype','path_in','path_out',...
    'nc_dim','nc_u','nc_v','nc_ssh','x_name','y_name','m_name','u_name','v_name','s_name',...
    'grid_ll','type_detection','extended_diags','streamlines','daystreamfunction','periodic',...
    'deg','resol','K','b','bx','Dx','Rd','gama','Rb','bi','bxi','Dxi','Rdi','gamai','Rb',...
    'stepF','T','f','g','dps','r','Dt','cut_off','D_stp','N_can','H','n_min',...
    'vel_epsil','k_vel_decay','dc_max','nR_lim','Np','nrho_lim')



