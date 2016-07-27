function [lon,lat,mask,u,v,ssh] = load_fields_AVISO(nc_dim,nc_u,nc_v,nc_ssh,b,res,deg)
%[lon,lat,mask,u,v,ssh]=load_field(nc_dim,nc_u,nc_v,nc_ssh,b,res {,deg});
%
%  Get netcdf and load the grid and the velocities field (also ssh if any)
%  for AVISO fields like. Other load_fields routines are used for models
%  or experiment.
%  if res>1 you need to call the funtion twice:
%   res=0: no interpolation.
%   res>1: interpolate the grid by a resolution factor res
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
%  For a description of the input parameters see param_eddy_tracking.m.
%
%-------------------------
%  Apr 2015 Briac Le Vu
%-------------------------
%
%=========================

if nargin==6
% No degradation by default
    deg = 1;
end

%----------------------------------------
% get netcdf
lon0 = double(ncread(nc_dim,'lon'))';
lat0 = double(ncread(nc_dim,'lat'))';
%mask0 = double(ncread(nc_dim,'mask'))'; ! mask0 not constant with time !
u0 = permute(ncread(nc_u,'u'),[2,1,3]);
v0 = permute(ncread(nc_v,'v'),[2,1,3]);

global type_detection
if type_detection>=2
    ssh0 = permute(ncread(nc_ssh,'ssh'),[2,1,3]);
else
    ssh0 = [];
end

%----------------------------------------
% produce degraded field 
lon0 = lon0(1:deg:end,1:deg:end);
lat0 = lat0(1:deg:end,1:deg:end);
u0 = u0(1:deg:end,1:deg:end,:);
v0 = v0(1:deg:end,1:deg:end,:);
if type_detection>=2
    ssh0 = ssh0(1:deg:end,1:deg:end,:);
end

% get the grid size
[N0,M0,stp] = size(u0);

%----------------------------------------
% ! AVISO mask changes with time !

% mask 2d
mask0 = min((u0(:,:,1:max(1,stp-7)).*v0(:,:,1:max(1,stp-7)))*0+1,[],3);
mask0(isnan(mask0)) = 0;

%lmask 3d
mask3d = repmat(mask0,[1,1,stp]);

%----------------------------------------
% fix fields to 0 in land and in miscellaneous AVISO mask
u0(mask3d==0) = 0;
v0(mask3d==0) = 0;
u0(isnan(u0)) = 0;
v0(isnan(v0)) = 0;
if type_detection>=2
    ssh0(mask3d==0) = NaN;
end

% Enlarge mask into land by b pixels

disp('"enlarge coastal mask" by adding b pixels of ocean to the coast')
mask1 = mask0;
for n=1:b
	mask2 = mask1;
    for i=1:N0
        for j=1:M0
            if mask1(i,j)==0
                if sum(sum(mask1(max(i-1,1):min(i+1,N0),max(j-1,1):min(j+1,M0))))~=0
                    mask2(i,j)=1;
if type_detection>=2
                    for k=1:stp
                        if n==1
                            ssh1 = ssh0(max(i-n,1):min(i+n,N0),max(j-n,1):min(j+n,M0),k);
                            ssh0(i,j,k) = nanmean(ssh1(:));
                        end
                    end
end
                end
            end
        end
    end
    mask1 = mask2;
end

%----------------------------------------
% Increase resolution r factor by linear interpolation

if res==1
	disp('NO INTERPOLATION')
    
    lon = lon0;
    lat = lat0;
    u = u0;
    v = v0;
    ssh = ssh0;
    mask = mask0;
    mask0 = mask1;
else
    disp(['"change resolution" by computing SPLINE INTERPOLATION res=',num2str(res)])

    % size for the interpolated grid
    N = fix(res*N0); % new size in lat
    M = fix(res*M0); % new size in lon

    % elemental spacing for the interpolated grid
    dlat = diff(lat0(1:2,1))/res;
    dlon = diff(lon0(1,1:2))/res;

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

    % initialize interp fields
    u = zeros([N M stp]);
    v = u;
    if type_detection>=2
        ssh = u;
    else
        ssh = [];
    end
    % Increase resolution of fields (interp2 with regular grid)
    for i=1:stp
        u(:,:,i) = interp2(lon0,lat0,squeeze(u0(:,:,i)),lon,lat,'*spline');
        v(:,:,i) = interp2(lon0,lat0,squeeze(v0(:,:,i)),lon,lat,'*spline');
        if type_detection>=2
            ssh(:,:,i) = interp2(lon0,lat0,squeeze(ssh0(:,:,i)),lon,lat,'linear');
        end
    end
end

disp(' ')

%----------------------------------------
% Mask velocities with the enlarged mask
mask3d = repmat(mask0,[1,1,stp]);
u(mask3d==0) = NaN;
v(mask3d==0) = NaN;

