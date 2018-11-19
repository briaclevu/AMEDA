function [x,y,mask,u,v,ssh] = load_fields_AVISO(stp,resolution,degra)
%[x,y,mask,u,v,ssh]=load_field(stp {,resolution,degra});
%
%  Get netcdf and load the grid and the velocities field (also ssh if any)
%  for AVISO fields like. Other load_fields routines are used for models
%  or experiment.
%  - degrad ('deg') is a sample factor (use in some experiment).
%  - resolution ('resol') is the factor to interpolate the grid by a resolution factor res.
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
%  For a description of the input parameters see param_eddy_tracking.m
%
%-------------------------
%  June 2016 Briac Le Vu
%-------------------------
%
%=========================

%----------------------------------------
% load keys_sources and parameters (use mod_eddy_params.m first)
load('param_eddy_tracking')

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
u0 = double(squeeze(permute(ncread(nc_u,u_name,[1 1 stp],[Inf Inf 1]),[2,1,3])));
v0 = double(squeeze(permute(ncread(nc_v,v_name,[1 1 stp],[Inf Inf 1]),[2,1,3])));
if type_detection>=2
    ssh0 = double(squeeze(permute(ncread(nc_ssh,s_name,[1 1 stp],[Inf Inf 1]),[2,1,3])));
else
    ssh = [];
end

%----------------------------------------
% produce degraded field 
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
if resol==1

    disp('NO INTERPOLATION')
    
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

else

    disp(['Change resolution by computing 2D SPLINE INTERPOLATION ',...
        'by a factor ',num2str(resol)])

    %----------------------------------------
    % Increase resolution of the grid
    % size for the interpolated grid
    Ni = resol*(N-1)+1; % new size in lat
    Mi = resol*(M-1)+1; % new size in lon
    % elemental spacing for the interpolated grid with regular grid
    dy = diff(y(1:2,1))/resol;
    dx = diff(x(1,1:2))/resol;
    % interpolated grid
    [xi,yi] = meshgrid((0:Mi-1)*dx+min(x(:)),(0:Ni-1)*dy+min(y(:)));

    %----------------------------------------
    % Increase resolution of the mask
    maski = interp2(x,y,mask,xi,yi);
    maski(isnan(maski) | maski < .5) = 0;
    maski(maski >= .5) = 1;
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
    % Increase resolution of the fields
    % fix fields to 0 in land
    u(mask==0 | isnan(u)) = 0;
    v(mask==0 | isnan(v)) = 0;
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
        sshi(maski1==0) = NaN;
    end

    %----------------------------------------
    % Export interpolated fields
    x = xi;
    y = yi;
    mask = maski;
    u = ui;
    v = vi;
    ssh = sshi;

end
