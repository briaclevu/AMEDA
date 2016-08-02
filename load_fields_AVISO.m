function [x,y,mask,u,v,ssh] = load_fields_AVISO(stp,resolution,degra)
%[x,y,mask,u,v,ssh]=load_field(stp {,resolution,degra});
%
%  Get netcdf and load the grid and the velocities field (also ssh if any)
%  for AVISO fields like. Other load_fields routines are used for models
%  or experiment.
%  degrad ('deg') is a sample factor (use in some experiment).
%  resolution ('resol') is the factor to interpolate the grid by a resolution factor res.
%
%  Enlarge the mask into land by adding b pixels of nil velocities into the
%  land and 1 pixel in land from an averaged of 9 neighbours for ssh
%	u(b)=0,u(~mask(b))=nan;
%	v(b)=0,v(~mask(b))=nan;
%	ssh(b)=mean(ssh)(b-1) or mean(ssh),ssh(~mask(b-1))=nan;
%
%  Output are matlab matrice used with tracking_plot routines and should be
%  saved in fields.mat file and fieldsi in fields_inter.mat.
%  fields size must be [lat,lon,time]
%  output velocities must be in m/s and ssh in m
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
x = lon0(1:deg:end,1:deg:end);
y = lat0(1:deg:end,1:deg:end);
u = u0(1:deg:end,1:deg:end);
v = v0(1:deg:end,1:deg:end);
if type_detection>=2
    ssh = ssh0(1:deg:end,1:deg:end);
end

% get the grid size
[N,M] = size(u);

% mask 2d ! AVISO mask changes with time ! use that one
%----------------------------------------
mask = u.*v*0+1;
mask(isnan(mask)) = 0;

% fix fields to 0 in land and in miscellaneous AVISO mask
%----------------------------------------
u(mask==0 | isnan(u)) = 0;
v(mask==0 | isnan(v)) = 0;

if type_detection>=2
    ssh(mask0==0) = NaN;
end

% Enlarge mask into land by b pixels and compute ssh in the first land pixel if needed
%----------------------------------------
mask1 = mask;
for n=1:max(b(:))
    mask2 = mask1;
    for i=1:N
        for j=1:M
            if mask1(i,j)==0 && max(max(b(max(i-1,1):min(i+1,N),max(j-1,1):min(j+1,M))))>=n &&...
                sum(sum(mask1(max(i-1,1):min(i+1,N),max(j-1,1):min(j+1,M))))~=0
                mask2(i,j)=1;
                if type_detection>=2 && n==1 && isnan(ssh(i,j))
                    ssh1 = ssh(max(i-n,1):min(i+n,N),max(j-n,1):min(j+n,M));
                    ssh(i,j) = nanmean(ssh1(:));
                end
            end
        end
    end
    mask1 = mask2;% enlarged mask
end

% Build mask ssh
%----------------------------------------
if type_detection>=2
    maskssh = ssh*0+1;
    maskssh(isnan(maskssh)) = 0;
end

% Mask velocities with the enlarged mask
%----------------------------------------
u(mask1==0) = NaN;
v(mask1==0) = NaN;

% Increase resolution r factor by linear interpolation
%----------------------------------------
if resol==1

    disp('NO INTERPOLATION')

else

    disp(['"change resolution" by computing 2D SPLINE INTERPOLATION res=',num2str(resol)])

    % size for the interpolated grid
    Ni = res*(N-1)+1; % new size in lat
    Mi = res*(M-1)+1; % new size in lon

    % elemental spacing for the interpolated grid
    dy = diff(y(1:2,1))/resol;
    dx = diff(x(1,1:2))/resol;

    % interpolated grid
    [xi,yi] = meshgrid([0:Mi-1]*dx+min(x(:)),[0:Ni-1]*dy+min(y(:)));

    % Increase resolution of the mask
    maski = interp2(x,y,mask,xi,yi);
    maski(isnan(maski) | maski < 1) = 0;
    
    if type_detection>=2
        % Increase resolution of the ssh mask
        maskissh = interp2(x,y,maskssh,xi,yi);
        maskissh(isnan(maskissh) | maskissh < 1) = 0;
    end

    % Compute the enlarged mask to the new resolution
    maski1 = interp2(x,y,mask1,xi,yi);
    maski1(isnan(maski1)) = 0;
    maski1(maski1 > 0) = 1;

    % Increase resolution of fields (interp2 with regular grid)
    ui = interp2(x,y,u,xi,yi,'*spline');
    vi = interp2(x,y,v,xi,yi,'*spline');
    if type_detection>=2
        ssh1 = get_missing_val_2d(x,y,ssh);
        sshi = interp2(x,y,ssh1,xi,yi,'*spline');
    else
        sshi = [];
    end
    
    % Mask velocities and ssh with their enlarged mask
    ui(maski1==0) = NaN;
    vi(maski1==0) = NaN;
    if type_detection>=2
        sshi(maskissh==0) = NaN;
    end

    % Export interpolated fields
    x = xi;
    y = yi;
    mask = maski;
    u = ui;
    v = vi;
    ssh = sshi;

end

