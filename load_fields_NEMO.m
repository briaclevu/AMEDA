function [lon,lat,mask,u,v,ssh] = load_fields_NEMO
%[lon,lat,mask,u,v,ssh]=load_field;
%
%  Get netcdf and load the grid and the velocities field (also ssh if any)
%  in case of res~=1: interpolate the grid by a resolution factor res
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

% load key_source and parameters
load('param_eddy_tracking')

% get netcdf
lon0 = double(ncread(nc_dim,lat_name))';
lat0 = double(ncread(nc_dim,lon_name))';
mask0 = double(ncread(nc_dim,mask_name))';
u0 = squeeze(permute(ncread(nc_u,u_name),[2,1,3,4]));
v0 = squeeze(permute(ncread(nc_v,v_name),[2,1,3,4]));

if type_detection>=2
    ssh0 = squeeze(permute(ncread(nc_ssh,ssh_name),[2,1,3,4]));
else
    ssh0 = [];
end

% get the grid size
[N0,M0] = size(lat0);
stp = size(u0,3);

%----------------------------------------
% mask 3d
mask3d = repmat(mask0,[1,1,stp]);

% fix fields to 0 in land
u0(mask3d==0) = 0;
v0(mask3d==0) = 0;
if type_detection>=2
    ssh0(mask3d==0) = NaN;
end

%----------------------------------------
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
    disp(['"change resolution" by computing LINEAR INTERPOLATION res=',num2str(res)])

    N = fix(res*N0); % new size in lat
    M = fix(res*M0); % new size in lon

    dlat = diff(lat0(1:2,1))/res;
    dlon = diff(lon0(1,1:2))/res;

    [lon,lat] = meshgrid(min(min(lon0)):dlon:(M-1)*dlon+min(min(lon0)),...
        min(min(lat0)):dlat:(N-1)*dlat+min(min(lat0)));

    % Increase resolution of the mask
    mask = griddata(lon0,lat0,mask0,lon,lat);
    mask(isnan(mask) | mask < 1) = 0;
    
    % Compute the enlarged mask to the new resolution
    mask0 = griddata(lon0,lat0,mask1,lon,lat);
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
    % Increase resolution of fields (griddata with irregular grid)
    for i=1:stp
        u(:,:,i) = griddata(lon0,lat0,squeeze(u0(:,:,i)),lon,lat);
        v(:,:,i) = griddata(lon0,lat0,squeeze(v0(:,:,i)),lon,lat);
        if type_detection>=2
            ssh(:,:,i) = griddata(lon0,lat0,squeeze(ssh0(:,:,i)),lon,lat);
        end
    end
end

disp(' ')

%----------------------------------------
% Mask velocities with the enlarged mask
mask3d = repmat(mask0,[1,1,stp]);
u(mask3d==0) = NaN;
v(mask3d==0) = NaN;
