function [x,y,mask,u,v,xi,yi,maski,ui,vi] = load_fields_ROMS(mat,b,box,deg,res)
%[x,y,mask,u,v,xi,yi,maski,ui,vi]=load_fields_inter(mat,b,box {,deg,res})
%
%  Load velocities field from a 2x2 kmxkm from a mat file
%
%  Extract degraded field by a factor 'deg', repeat 'box' ipxels if periodic
%  add sponge layer with 'b' pixels
%  and then interpolate by a factor 'res'
%
%  In case of res~=1: interpolate the grid by a resolution factor res
%
%  Output are matlab matrice used with eddy_tracking routines
%
%  Fields size must be [y,x,time]
%
%  Output velocities must be in m/s
%
%  For a description of the input parameters see param_eddy_tracking.m
%
%-------------------------
%  Apr 2015 Briac Le Vu
%-------------------------
%
%=========================

global periodic

if nargin==4
% No interpolation by default
    res = 1;
elseif nargin==3
% No degradation and No interpolation by default
    deg = 1;
    res = 1;
end

%% load fields
load(mat)
uroms = double(permute(uroms,[2,1,3]));
vroms = double(permute(vroms,[2,1,3]));

% get the u,v grids sizes
[Nu,Mu,P] = size(uroms);
[Nv,Mv,~] = size(vroms);

% regular initial T grid and mask
maskr = ones(Nv,Mu);
[xr,yr] = meshgrid((1/2:Mu-1/2)*xd,(1/2:Nv-1/2)*xd);

%% Enlarge degraded field
if periodic
    % add box pixels to the boundaries (East-West periodicity)
    % u and v matrices are expanded by adding the ending columns
    % to the beginning of the domain, and the beginning columns at its end.
    borders = box;
    % enlaregd uroms with the same values at West and East boundary
    uenl=[uroms(:,end-box:end-1,:),uroms,uroms(:,2:box+1,:)];
    % enlarged vroms with the same 2 first values at the West and the East
    venl=[vroms(:,end-box-1:end-2,:),vroms,vroms(:,3:box+2,:)];
else
    % add no pixels to the West and East boundaries
    borders = 0;
    uenl=uroms;
    venl=vroms;
end

% u enlarged grid
[xu,yu,zu] = meshgrid((1/2-borders:Mu-1/2+borders)*xd,(0:Nu-1)*xd,1:P);
% v enlarged grid
[xv,yv,zv] = meshgrid((-borders:Mv-1+borders)*xd,(1/2:Nv-1/2)*xd,1:P);

% common enlarged degraded grid interpolated in the middle of y and x grid
[xq,yq,zq] = meshgrid((1/2-borders:deg:Mu-1/2+borders)*xd,(1/2:deg:Nv-1/2)*xd,1:P);

% degraded enlarged degraded field interpolate u and v on the middle
uq = interp3(xu,yu,zu,uenl,xq,yq,zq);
vq = interp3(xv,yv,zv,venl,xq,yq,zq);

% get the enlarged degraded grid size
[Nq,Mq,~] = size(xq);

% produce common enlarged mask
maskq = ones(Nq,Mq);

%% Add a sponge layer for LNAM computation to the enlarged degraded fields
% add b pixel to the East and West boundaries with value 0
uq = [zeros(Nq,b,P),uq,zeros(Nq,b,P)];
vq = [zeros(Nq,b,P),vq,zeros(Nq,b,P)];

% add b pixels to the North and South boundaries with value 0
u = [zeros(b,Mq+2*b,P);uq;zeros(b,Mq+2*b,P)];
v = [zeros(b,Mq+2*b,P);vq;zeros(b,Mq+2*b,P)];

% produce x and y enlarged degraded grid with sponge (2x2 kmxkm)
[x,y] = meshgrid(((1/2-borders-b*deg):deg:(Mu-1/2+borders+b*deg))*xd,((1/2-b*deg):deg:(Nv-1/2+b*deg))*xd);

% size degraded field
[N,M] = size(x);

% produce mask of the initial grid on enlarged degraded grid
% !!! probleme during compute_psi which use the mask !!!
%mask = interp2(xr,yr,maskr,x,y);
%mask(isnan(mask) | mask < 1) = 0;

% mask of enlarged grid
mask = interp2(xq(:,:,1),yq(:,:,1),maskq,x,y);
mask(isnan(mask) | mask < 1) = 0;

%% interpolation on interpolated grid by 'res'
if res==1
	disp('NO INTERPOLATION')
    
    xi = x;
    yi = y;
    ui = u;
    vi = v;
    maski = mask;
    
else
    disp(['"change resolution" by computing SPLINE INTERPOLATION res=',num2str(res)])
    
    % size for the interpolated grid
    Ni = fix(res*N); % new size in y
    Mi = fix(res*M); % new size in x
    
    % elemental spacing for the interpolated grid
    dx = diff(x(1,1:2))/res;
    dy = diff(y(1:2,1))/res;

    % interpolated grid
    [xi,yi] = meshgrid(min(x(:)):dx:(Mi-1)*dx+min(x(:)),...
            min(y(:)):dy:(Ni-1)*dy+min(y(:)));

    % Increase resolution of the mask
    maski = interp2(x,y,mask,xi,yi);
    maski(isnan(maski) | maski < 1) = 0;
    
    % initialize interp fields
    ui = zeros([Ni Mi P]);
    vi = ui;

    % Increase resolution of fields (interp2 with regular grid)
    for i=1:P
        ui(:,:,i) = interp2(x,y,squeeze(u(:,:,i)),xi,yi,'spline');
        vi(:,:,i) = interp2(x,y,squeeze(v(:,:,i)),xi,yi,'spline');
    end

end

disp(' ')

