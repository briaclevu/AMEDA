function [lon,lat,mask,u,v,ssh] = load_fields_AVISO(stp,resolution,degra)
%[lon,lat,mask,u,v,ssh]=load_field(stp {,resolution,degra});
%
%  Get netcdf and load the grid and the velocities field (also ssh if any)
%  for AVISO fields like.
%  degrad ('deg') is a sample factor (use in some experiment).
%  resolution ('resol') is the factor to interpolate the grid (to find centers)
%
%  Enlarge the mask into land by adding b pixels of nil velocities into the
%  land and 1 pixel in land from an averaged of 9 neighbours for ssh
%	u(b)=0,u(~mask(b))=nan;
%	v(b)=0,v(~mask(b))=nan;
%	ssh(b)=mean(ssh)(b-1) or mean(ssh),ssh(~mask(b-1))=nan;
%
%  Output are matlab matrice used with tracking_plot routines and should be
%  saved in fields.mat file. fields size must be [lat,lon,time]
%  output velocities must be in m/s and ssh in m
%
%  Other load_fields routines are used for models or experiment.
%
%  For a description of the input parameters see param_eddy_tracking.m.
%
%-------------------------
%  June 2016 Briac Le Vu
%-------------------------
%
%=========================

% load keys_sources and parameters (use mod_eddy_params.m first)
%----------------------------------------
load('param_eddy_tracking')

% replace parameters by arguments
%----------------------------------------
if nargin==3
    deg = degra;
    resol = resolution;
elseif nargin==2
    resol = resolution;
end

% get netcdf
%----------------------------------------
lon0 = double(ncread(nc_dim,x_name))';
lat0 = double(ncread(nc_dim,y_name))';
% ! mask0 not constant with time ! don't trust that one
%mask0 = double(ncread(nc_dim,m_name))';
u0 = squeeze(permute(ncread(nc_u,u_name,[1 1 stp],[Inf Inf 1]),[2,1,3]));
v0 = squeeze(permute(ncread(nc_v,v_name,[1 1 stp],[Inf Inf 1]),[2,1,3]));

if type_detection>=2
    ssh0 = squeeze(permute(ncread(nc_ssh,s_name,[1 1 stp],[Inf Inf 1]),[2,1,3]));
else
    ssh0 = [];
end

% produce degraded field 
%----------------------------------------
lon0 = lon0(1:deg:end,1:deg:end);
lat0 = lat0(1:deg:end,1:deg:end);
u0 = u0(1:deg:end,1:deg:end);
v0 = v0(1:deg:end,1:deg:end);
if type_detection>=2
    ssh0 = ssh0(1:deg:end,1:deg:end);
end

% get the grid size
[N0,M0] = size(u0);

% mask 2d ! AVISO mask changes with time ! use that one
%----------------------------------------
mask0 = u0.*v0*0+1;
mask0(isnan(mask0)) = 0;

% fix fields to 0 in land and in miscellaneous AVISO mask
%----------------------------------------
u0(mask0==0) = 0;
v0(mask0==0) = 0;
u0(isnan(u0)) = 0;
v0(isnan(v0)) = 0;
if type_detection>=2
    ssh0(mask0==0) = NaN;
end

% Enlarge mask into land by b pixels
%----------------------------------------
mask1 = mask0;
for n=1:b
	mask2 = mask1;
    for i=1:N0
        for j=1:M0
            if mask1(i,j)==0
                if sum(sum(mask1(max(i-1,1):min(i+1,N0),max(j-1,1):min(j+1,M0))))~=0
                    mask2(i,j)=1;
                    if type_detection>=2 && n==1
                        ssh1 = ssh0(max(i-n,1):min(i+n,N0),max(j-n,1):min(j+n,M0));
                        ssh0(i,j) = nanmean(ssh1(:));
                    end
                end
            end
        end
    end
    mask1 = mask2;
end

% Increase resolution r factor by linear interpolation
%----------------------------------------
if resol==1
    lon = lon0;
    lat = lat0;
    u = u0;
    v = v0;
    ssh = ssh0;
    mask = mask0;
    mask0 = mask1;
else
    % size for the interpolated grid
    N = fix(resol*N0); % new size in lat
    M = fix(resol*M0); % new size in lon

    % elemental spacing for the interpolated grid
    dlat = diff(lat0(1:2,1))/resol;
    dlon = diff(lon0(1,1:2))/resol;

    % interpolated grid
    [lon,lat] = meshgrid(min(lon0(:)):dlon:(M-1)*dlon+min(lon0(:)),...
        min(lat0(:)):dlat:(N-1)*dlat+min(lat0(:)));

    % Increase resolution of the mask
    mask = interp2(lon0,lat0,mask0,lon,lat);
    mask(isnan(mask) | mask < 1) = 0;
    
    % Compute the enlarged mask to the new resolution
    mask0 = interp2(lon0,lat0,mask1,lon,lat);
    mask0(isnan(mask0)) = 0;
    mask0(mask0 > 0) = 1;

    % Increase resolution of fields (interp2 with regular grid)
    u = interp2(lon0,lat0,u0,lon,lat,'*spline');
    v = interp2(lon0,lat0,v0,lon,lat,'*spline');
    if type_detection>=2
        ssh = interp2(lon0,lat0,ssh0,lon,lat,'linear');
    else
        ssh = [];
    end
end

disp(' ')

% Mask velocities with the enlarged mask
%----------------------------------------
u(mask0==0) = NaN;
v(mask0==0) = NaN;

