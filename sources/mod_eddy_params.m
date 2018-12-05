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
%   - Rd: First Baroclinic Rossby Radius of Deformation use as the typical
%       size of serached eddies
%   - gama: resolution coefficient which is the number of pixels per Rd.
%   - resol: factor of interpolation of the fields (velocity and ssh)
%       to compute center detection. 2 or 3 seems reasonable in term of
%       computation time for 1/8° AVISO fields.
%   - f: The Coriolis parameter
%   - x and y for the native grid and xi, yi for the regular grid
%
% Fixed parameters
% (you can play with Dt and cut_off but don't touch too much others):
%   - K: LNAM(LOW<0) threshold to detect the potential eddy centers
%   - V_eddy: center maximale speed during tracking of an eddy
%   - Dt: delay in days for the tracking tolerance
%   - cut_off: in days for eddies duration filtration after the tracking
%           0 : use the turnover time from each eddies
%           1 : keep all the tracked eddies
%   - D_stp: in steps to define the averaged eddy features during the tracking
%   - N_can: maximum number of candidat the assignment must resolve
%   - DH: ssh space between scanned streamlines in the initial area
%       for each eddies (see 'box' parameter below)
%   - n_min: minimal number of point to defined a contour as streamlines
%   - vel_epsil: minimal difference between 2 tests of the velocity
%       to admit an increase
%   - k_vel_decay: coefficient of velmax to detect a decrease
%   - nR_lim: limite of the eddy size in temr of Rd
%   - Np: minimal number of point to calculate curvature along a segment
%   - nrho_lim: limite of the negative curvature for a single eddy
%       interaction (merging or splitting)
%   - dc_max: maximal distance between 2 eddies contours
%   - lat_min: minimal latitude to scan potential centers among LNAM extrema 
%
%----------------------------------------------
% This parametrisation come from the test with AVISO, ROMS and PIV
% see LE VU et al. 2018
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
if strcmp(source,'NEMO')
    mask0 = squeeze(double(ncread(nc_dim,m_name,[1 1 level 1],[Inf Inf 1 1])))';
elseif strcmp(source,'py')
    mask0 = ones(size(lat0))';
else
    mask0 = squeeze(double(ncread(nc_dim,m_name,[1 1],[Inf Inf])))';
end

% produce degraded field 
if strcmp(source,'py')
    x = lon0(1:deg:end,1:deg:end)/1000;
    y = lat0(1:deg:end,1:deg:end)/1000;
else
    x = lon0(1:deg:end,1:deg:end);
    y = lat0(1:deg:end,1:deg:end);
end
mask = mask0(1:deg:end,1:deg:end);

[N,M] = size(x);

% Meshgrid size
if grid_ll
    if grid_reg
        Dx = get_Dx_from_ll(x,y);
    else
        disp('Irregular grid : regridding')
        dy = mean( mean( diff(y,[],1) ));
        dx = mean( mean( diff(x,[],2) ));
        % regridded grid
        [xr, yr] = meshgrid(min(x(:)):dx:max(x(:)+dx),min(y(:)):dy:max(y(:))+dy);
        [N,M] = size(xr);
        Dx = get_Dx_from_ll(xr,yr);
        % regridded mask
        maskr = griddata(x,y,mask,xr,yr);
        maskr(isnan(maskr) | maskr < .5) = 0;
        maskr(maskr >= .5) = 1;
        maskr(xr(:)<min(x(:)) | xr(:)>max(x(:)) | yr(:)<min(y(:)) | yr(:)>max(y(:))) = 0;
    end
else
    Dx = ( abs( [diff(x,1,2) x(:,end)-x(:,end-1)] ) + abs( [diff(y,1,1);y(end,:)-y(end-1,:)] ) )/2;
end

% Calculate coriolis parameter
if grid_ll
    if grid_reg
        f = 4*pi/T*sind(y); % in s-1
    else
        f = 4*pi/T*sind(yr); % in s-1
    end
else
    f = 4*pi/T*(y*0+1);
end

g = 9.8; % m.s-2

%% fixed parameters
%----------------------------------------------

% Scanning parameters:
%----------------------------------------------
% H is the ssh space of streamlines between min and max ssh in the box
% to be scanned during detection of potential centers (centers0 to centers)
% and shapes (centers to centers2).
% After test [0.001 - 0.002] seems a good compromise.
% This parameter is directly link to the time computation!
DH = 0.002; % (m) ssh space

% H_lim is the maximal number of streamlines scanned in a domain around a
% center
nH_lim = 200;

% these contours must contain a minimal number of dots ([4-6])
n_min = 6; % minimal number of point to defined a contour as streamlines

% minimal difference increase ([1-2]%) between 2 tests of the velocity
% or deta
epsil = 1e-2; % %

% coefficient of velmax to detect a decrease ([0.95-0.99])
k_vel_decay = 0.97; % in term of velmax

% size limite in term of deformation radius ([4-7]Rd) for a single eddy
nR_lim = 100; % [4-7]Rd
% !!! decide to no limit the size which introduce a bias
% by artificially limiting size of gyre!!!

% minimal number of point to calculate curvature along a segment at point i
Np = 3; % (segment i-Np:i+Np)

% limite of the length of the contour with negative curvature for a single
% eddy from empirical determination based on ROMS, PIV and AVISO tests
nrho_lim = 0.2; % [1/5-1/2]*2*pi*Rmax of the equivalent circle

% minimal latitude used in terrestrial coordiante (grid_ll=1) needed
% to search for potential centers amon the LNAM extrema in mod_eddy_centers
lat_min = 5; % [1-15]

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
% typical eddy speed deplacement in a day
% 0 to use <Vmax> for each eddy instead
V_eddy = 6.5; % km/day

% maximal delay after its last detection for tracking an eddy [1-15]
% represent the temporal correlation depending on the coverage of an area
% 2 steps for model, 3-5 steps for imagery experiment.
%!!! in case of AVISO this will be ajusted with the error map of aviso !!!
if strcmp(source,'AVISO')
    Dt = 10; % in days
else
    Dt = 2; % days
end

% minimal duration of recorded eddies [0-100]
% 0 for keeping only eddies longer than twice their turnover time
cut_off = 0; % in days

% number of past steps to define the averaged eddy features during the tracking
if strcmp(source,'AVISO')
    D_stp = 4; % in steps
else
    D_stp = 2; % in steps
end

% number of candidats during the tracking (better to be high)
N_can = 30;

%% Calculate parameters needed for AMEDA at Dx and a given Rd
%----------------------------------------------

% Interpolate First Baroclinic Rossby Radius of Deformation
% (taken from 2D file Rossby_radius computed using Chelton et al. 1998)
%----------------------------------------------
if exist([mat_Rd,'.mat'],'file')
    load(mat_Rd)
    eval([name_Rd,'(',name_Rd,'<',num2str(Rd_typ/3),')=',num2str(Rd_typ/3),';']);
    if grid_reg
        eval(['Rd0 = griddata(lon_Rd,lat_Rd,',name_Rd,',x,y);']) % 10km in average AVISO 1/8
        Rd = Rd0;
        Rd(mask==0)=nan;
    else
        eval(['Rd0 = griddata(lon_Rd,lat_Rd,',name_Rd,',xr,yr);']) % 10km in average AVISO 1/8
        Rd = Rd0;
        Rd(maskr==0)=nan;
    end
end

if ~exist([mat_Rd,'.mat'],'file') || isnan(nanmean(Rd(:))) || ~exist('Rd','var')
    disp('!!! WARNING or ERROR !!! Rossby Radius computation impossible - check the Rd input file')
    disp(['AMEDA will use typical a constant typical radius of ',num2str(Rd_typ),' km instead'])
    if grid_reg
        Rd = x * 0 + Rd_typ;
    else
        Rd = xr * 0 + Rd_typ;
    end
end

% Resolution parameters:
%----------------------------------------------
% gama is resolution coefficient which is the number of pixels per Rd.
% After test gama>2.4 is required to get the max number of eddies.
gama = Rd ./ Dx; % [0.1-1.5] for AVISO 1/8 (0.8 in average)

% resol is an integer used to improve the precision of centers detection
% close to 3 pixels per Rd. resol can goes up to 3
resol = max(1,min(3,round(3/nanmean(gama(:))))); % [1 - 3]

% Detection parameters (cf. Le Vu et al. 2018 for test results):
%----------------------------------------------
% K is LNAM(LOW<0) threshold to fixed contours where to detect potential center
% (one per contour). After test: no sensibility to K between 0.2 and 0.7
% but 0.7 give better time performance
K = 0.7; % [.4-.8]

% b is half length of the box in pixels used to normalise the LNAM and LOW.
% After test the optimal length of the box ( Lb = 2b*Dx*deg ) is fixed
% to 1.2 the size of Rd (Lb/Rd=1.2) to detect a maximum number of eddies.
b = max(1,round((1.2*gama)/2)); % always 1 for AVISO 1/8

% Rb (=Lb/Rd) is to check that the b setting is optimal for a given gama.
% the optimal Rb stay close to 1 (0.5 - 1.5).
% !!! Because optimal b is directly linked to gama, you start missing
% smaller eddies when gama is far from 2.4 at b=1 or gama far from 5 at b=2.
% (e.g. for AVISO 1/8 (gama~0.8) we have Rb~2.5) !!!
Rb = 2*b ./ gama; % [1-100] for AVISO 1/8 (2.5 in average)

% box is half length of the box used in pixels to find contour streamlines
% in pixels box = 5 times the Rd is enough to start testing any eddies
bx = max(1,round(2*gama*5)); % [1-14] for AVISO 1/8 (8 in average)

%% interpolated parameters
%----------------------------------------------

if resol==1

    if grid_reg
        xi = x;
        yi = y;
        maski = mask;
    else
        xi = xr;
        yi = yr;
        maski = maskr;
    end
    f_i = f;
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
    if grid_reg
        dy = diff(y(1:2,1))/resol;
        dx = diff(x(1,1:2))/resol;
    else
        dy = diff(yr(1:2,1))/resol;
        dx = diff(xr(1,1:2))/resol;
    end

    % interpolated grid
    [xi,yi] = meshgrid((0:Mi-1)*dx+min(x(:)),(0:Ni-1)*dy+min(y(:)));

    % regridded mask
    if grid_reg
        maski = griddata(x,y,mask,xi,yi);
    else
        maski = griddata(xr,yr,maskr,xi,yi);
        maski(xi(:)<min(x(:)) | xi(:)>max(x(:)) | yi(:)<min(y(:)) | yi(:)>max(y(:))) = 0;
    end
    maski(isnan(maski) | maski < .5) = 0;
    maski(maski >= .5) = 1;
    
    % Compute interpolated b and bx
    if grid_reg
        f_i = griddata(x,y,f,xi,yi)*resol;
        bi = round(griddata(x,y,b,xi,yi))*resol;
        bxi = round(griddata(x,y,bx,xi,yi))*resol;
    else
        f_i = griddata(xr,yr,f,xi,yi)*resol;
        bi = round(griddata(xr,yr,b,xi,yi))*resol;
        bxi = round(griddata(xr,yr,bx,xi,yi))*resol;
    end

    % Dx, Rd, gama on interpolated grid
    if grid_reg
        Rdi = griddata(x,y,Rd,xi,yi); % 10km in average AVISO 1/8
        Dxi = griddata(x,y,Dx,xi,yi);
    else
        Rdi = griddata(xr,yr,Rd,xi,yi); % 10km in average AVISO 1/8
        Dxi = griddata(xr,yr,Dx,xi,yi);
    end
    gamai = Rdi ./ Dxi;

end
    
%% Save parameters and paths
%----------------------------------------------
save([path_out,'param_eddy_tracking'],...
    'postname','config','runname','sshtype','path_in','path_out','path_rossby','level',...
    'nc_dim','nc_u','nc_v','nc_ssh','x_name','y_name','m_name','u_name','v_name','s_name',...
    'grid_ll','grid_reg','type_detection','extended_diags','streamlines','daystreamfunction','periodic','nrt',...
    'deg','resol','level','Rd_typ','K','b','bx','Dx','Rd','gama','Rb','bi','bxi','Dxi','Rdi','gamai','Rb',...
    'stepF','T','f','f_i','g','dps','V_eddy','Dt','cut_off','D_stp','N_can','DH','nH_lim','n_min','lat_min',...
    'epsil','k_vel_decay','dc_max','nRmin','nR_lim','Np','nrho_lim')

% save native grid on interpolated regular grid
if ~grid_reg
    save([path_out,'gridvel_deg',num2str(deg),'_resol',num2str(resol)], 'x', 'y','mask','xi', 'yi','maski')
end


