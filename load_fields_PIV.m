function [x,y,mask,u,v,xi,yi,maski,ui,vi] = load_fields_PIV(mat,b,deg,res)
%[x,y,mask,u,v,xi,yi,maski,ui,vi]=load_fields_inter(mat,b {,deg,res})
%
%  Load velocities field from a 2.18x2.20 mm from mat file
%
%  Extract degraded field by a factor 'deg',
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
%  Nov 2015 Briac Le Vu
%-------------------------
%
%=========================

if nargin==3
% No interpolation by default
    res = 1;
elseif nargin==2
% No degradation and No interpolation by default
    deg = 1;
    res = 1;
end

%% load fields
load(mat)
xp = x;
yp = y;
maskp = mask;
up = double(permute(u,[1,2,3]));
vp = double(permute(v,[1,2,3]));

%% Degraded field
% produce degraded x and y around 0°N 
xq = xp(1:deg:end,1:deg:end);
yq = yp(1:deg:end,1:deg:end);
maskq = maskp(1:deg:end,1:deg:end);

% extract degraded field u and v
uq = up(1:deg:end,1:deg:end,:);
vq = vp(1:deg:end,1:deg:end,:);

% get the degraded grid size
[Nq,Mq,P] = size(uq);

%% Add a sponge layer for LNAM computation to the enlarged degraded fields
% add b pixel to all boundaries 
[x,y] = meshgrid([xq(1,1)-(b:-1:1)*dx*deg,xq(1,:),xq(1,end)+(1:b)*dx*deg],...
	[yq(1,1)-(b:-1:1)*dy*deg,yq(:,1)',yq(end,1)+(1:b)*dy*deg]);

% add b pixels to the bundaries (0 east-west)
maskq = [zeros(Nq,b),maskq,zeros(Nq,b)];
vq = [zeros(Nq,b,P),vq,zeros(Nq,b,P)];
uq = [zeros(Nq,b,P),uq,zeros(Nq,b,P)];

% add b pixels to the bundaries (0 north-south)
mask = [zeros(b,Mq+2*b);maskq;zeros(b,Mq+2*b)];
v = [zeros(b,Mq+2*b,P);vq;zeros(b,Mq+2*b,P)];
u = [zeros(b,Mq+2*b,P);uq;zeros(b,Mq+2*b,P)];

% size degraded field
[N,M] = size(x);

%% Interpolation on interpolated grid by 'res'
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

