function [x,y,mask,u,v,ssh] = load_fields_CROCO(stp, resolution, degra)
%[x,y,mask,u,v,ssh]=load_field(stp {,resolution,degra});
%
%  Get netcdf and load the grid and the velocities field (also ssh if any)
%  - degrad ('deg') is a sample factor (use in some experiment).
%  - resolution ('resol') is the factor to interpolate the grid by a resolution factor res.
%   res=0: interpolation on regular grid if necesseray
%   res>1: interpolate the grid by a resolution factor res
%
%  Enlarge the mask into land by adding b pixels of nil velocities into the
%  land and 1 pixel in land from an averaged of 9 neighbours for ssh
%	u(b)=0,u(~mask(b))=nan;
%	v(b)=0,v(~mask(b))=nan;
%	ssh(b)=mean(ssh)(b-1) or mean(ssh),ssh(~mask(b-1))=nan;
%
%  Output are matlab matrice used with tracking_plot routines and should be
%  saved in fields.mat file and in fields_inter.mat.
%  fields size must be [lat,lon,time]
%  output velocities must be in m/s and ssh in m
%
%  For a description of the input parameters see mod_eddy_param.m.
%
%-------------------------
%  June 2018 Briac Le Vu and Romain Pennel
%-------------------------
%
%=========================

%----------------------------------------
% load keys_sources and parameters (use mod_eddy_params.m first)
load('param_eddy_tracking')
% critical depth during the simulation
hc=double(ncread(nc_dim,'hc'));
theta_s=double(ncread(nc_dim,'theta_s'));
theta_b=double(ncread(nc_dim,'theta_b'));
%hc = 350;theta_s=0;theta_b=0;% default MEDIONE

%----------------------------------------
% replace parameters by arguments
if nargin==3
% No interpolation by default
    deg = degra;
    resol = resolution;
elseif nargin==2
% No degradation and No interpolation by default
    resol = resolution;
end

%----------------------------------------
% get netcdf
disp(['Get grid and velocities field at step ',num2str(stp),' ...'])

lon0 = double(ncread(nc_dim,x_name))';
lat0 = double(ncread(nc_dim,y_name))';
mask0 = double(ncread(nc_dim,m_name))';
z0 = get_z_croco(nc_u,nc_ssh,stp,theta_s,theta_b,hc,0,'r',[1,3,2]);

%----------------------------------------
% interpolate u and v
u3d = double(squeeze(permute(ncread(nc_u,u_name,[1,1,1,stp],[Inf,Inf,Inf,1]),[4,3,2,1])));
v3d = double(squeeze(permute(ncread(nc_v,v_name,[1,1,1,stp],[Inf,Inf,Inf,1]),[4,3,2,1])));
lon_u = double(ncread(nc_dim,'nav_lon_u'))';
lat_u = double(ncread(nc_dim,'nav_lat_u'))';
lon_v = double(ncread(nc_dim,'nav_lon_v'))';
lat_v = double(ncread(nc_dim,'nav_lat_v'))';
zu = (z0(:,:,1:end-1)+z0(:,:,2:end))/2;
zv = (z0(:,1:end-1,:)+z0(:,2:end,:))/2;
if type_detection>=2
    ssh0 = squeeze(permute(ncread(nc_ssh,s_name,[1,1,stp],[Inf,Inf,1]),[2,1,3,4]));
else
    ssh = [];
end

%----------------------------------------
% get 2d fields at 'level' depth
disp(['Get velocities field  at the level ',num2str(level),'m at step ',num2str(stp),' ...'])
ub = vinterp(u3d,zu,level);
vb = vinterp(v3d,zv,level);
ub(isnan(ub))=0;
vb(isnan(vb))=0;

%----------------------------------------
% put velocities on rho grid
u0 = interp2(lon_u,lat_u,ub,lon0,lat0,'*spline');
v0 = interp2(lon_v,lat_v,vb,lon0,lat0,'*spline');

clear zu zv ub vb lon_u lat_u lon_v lat_v u3d v3d

%----------------------------------------
% work on degraded field 
if deg~=1
    disp(['  Fields are degraded by a factor ',num2str(deg)])
    disp('  (degraded grid becomes native grid)')
end

x = lon0(1:deg:end,1:deg:end);
y = lat0(1:deg:end,1:deg:end);
mask = mask0(1:deg:end,1:deg:end);
u = u0(1:deg:end,1:deg:end);
v = v0(1:deg:end,1:deg:end);
if type_detection>=2
    ssh = ssh0(1:deg:end,1:deg:end);
end

% get the grid size
[N,M] = size(x);

%----------------------------------------
% Increase resolution r factor by linear interpolation
if resol==1 && grid_reg

    disp('NO INTERPOLATION')
    
    % fix fields to NaN in land 
    u(mask==0) = NaN;
    v(mask==0) = NaN;

    % Enlarge mask into land by 1 pixel and compute ssh in the first land pixel if needed
    disp('Enlarge coastal mask by adding 1 pixel of ocean to the coast ...')
    for i=1:N
        for j=1:M
            if mask(i,j)==0 &&...
                    sum(sum(mask(max(i-1,1):min(i+1,N),max(j-1,1):min(j+1,M))))~=0
                u(i,j) = 0;
                v(i,j) = 0;
                if type_detection>=2 && isnan(ssh(i,j))
                    ssh1 = ssh(max(i-1,1):min(i+1,N),max(j-1,1):min(j+1,M));
                    ssh(i,j) = nanmean(ssh1(:));
                end
            end
        end
    end
    
elseif ( resol==1 && ~grid_reg ) || resol ~= 1

    if grid_reg
        disp(['Change resolution by computing 2D SPLINE INTERPOLATION ',...
        'by a factor ',num2str(resol)])
    else
        disp('No change in resolution, REGRIDDING from Arakawa to regular grid')
    end
    
    %----------------------------------------
    % Increase resolution of the mask
    [Ni,Mi] = size(xi);
    if size(x,1) ~= size(u,1) && Ni ~= resol*(N-1)+1
        disp(['Need to adapt gridvel.mat to deg=',num2str(deg),' and res=',num2str(resol)])
        return
    end
    
    %----------------------------------------
    % Enlarge mask into land by 1 pixel and compute ssh in the first land pixel if needed
    maski1 = maski;% enlarged mask
    disp('Enlarge coastal mask by adding 1 pixel of ocean to the coast ...')
    for i=1:Ni
        for j=1:Mi
            if maski(i,j)==0 &&...
                    sum(sum(maski(max(i-1,1):min(i+1,Ni),max(j-1,1):min(j+1,Mi))))~=0
                maski1(i,j)=1;
            end
        end
    end

    %----------------------------------------
    % fix fields to 0 in land 
    u(mask==0 | isnan(u)) = 0;
    v(mask==0 | isnan(v)) = 0;
    % Increase resolution of fields (griddata with irregular grid)
    ui = griddata(x,y,u,xi,yi,'cubic');
    vi = griddata(x,y,v,xi,yi,'cubic');
    if type_detection>=2
        sshi = griddata(x,y,ssh,xi,yi,'cubic');
    else
        sshi = [];
    end
    % Mask velocities and ssh with their enlarged mask
    ui(maski1==0) = NaN;
    vi(maski1==0) = NaN;
    if type_detection>=2
        sshi(maski1==0) = NaN;
    end

    %----------------------------------------
    % Export interpolated fields
    x = xi;
    y = yi;
    u = ui;
    v = vi;
    mask = maski;
    ssh = sshi;
       
end
