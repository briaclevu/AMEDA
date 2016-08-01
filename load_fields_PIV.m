<<<<<<< HEAD
function [x,y,mask,u,v,xi,yi,maski,ui,vi] = load_fields_PIV(mat,b0,bx0,res,deg)
%[x,y,mask,u,v,xi,yi,maski,ui,vi]=load_fields_inter(mat,b,bx {,res,deg})
=======
function [x,y,mask,u,v,xi,yi,maski,ui,vi] = load_fields_PIV(mat,b,deg,res)
%[x,y,mask,u,v,xi,yi,maski,ui,vi]=load_fields_inter(mat,b {,deg,res})
>>>>>>> ameda_v2
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

<<<<<<< HEAD
global path_out
global runname

if nargin==4
% No degradation by default
    deg = 1;
elseif nargin==3
=======
if nargin==3
% No interpolation by default
    res = 1;
elseif nargin==2
>>>>>>> ameda_v2
% No degradation and No interpolation by default
    deg = 1;
    res = 1;
end

%% load fields
load(mat)
<<<<<<< HEAD
xp = x';
yp = y';
maskp = mask';
up = double(permute(u,[2,1,3]));
vp = double(permute(v,[2,1,3]));

%% Degraded field
% produce degraded x and yas North hemispher rotation
%xq = xp(end:-deg:1,1:deg:end);
%yq = -yp(end:-deg:1,1:deg:end);
%maskq = maskp(end:-deg:1,1:deg:end);
=======
xp = x;
yp = y;
maskp = mask;
up = double(permute(u,[1,2,3]));
vp = double(permute(v,[1,2,3]));

%% Degraded field
% produce degraded x and y around 0°N 
>>>>>>> ameda_v2
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
<<<<<<< HEAD
[x,y] = meshgrid([xq(1,1)-(b0:-1:1)*dx*deg,xq(1,:),xq(1,end)+(1:b0)*dx*deg],...
    [yq(1,1)-(b0:-1:1)*dy*deg,yq(:,1)',yq(end,1)+(1:b0)*dy*deg]);

% add b pixels to the bundaries (0 east-west)
maskq = [zeros(Nq,b0),maskq,zeros(Nq,b0)];
vq = [zeros(Nq,b0,P),vq,zeros(Nq,b0,P)];
uq = [zeros(Nq,b0,P),uq,zeros(Nq,b0,P)];

% add b pixels to the bundaries (0 north-south)
mask = [zeros(b0,Mq+2*b0);maskq;zeros(b0,Mq+2*b0)];
v = [zeros(b0,Mq+2*b0,P);vq;zeros(b0,Mq+2*b0,P)];
u = [zeros(b0,Mq+2*b0,P);uq;zeros(b0,Mq+2*b0,P)];
=======
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
>>>>>>> ameda_v2

% size degraded field
[N,M] = size(x);

<<<<<<< HEAD
% 2D parameters
b    = x*0 + b0;
bx   = x*0 + bx0;
Rd   = x*0 + Rd(1,1);
gama = Rd ./ (Dx*deg);

%% Interpolation on interpolated grid by 'res'
if res==1

    disp('NO INTERPOLATION')
    
    xi = x;
    yi = y;
    maski = mask;
    ui = u;
    vi = v;
    bi = b;
    bxi = bx;
    Dxi = Dx;
    Rdi = Rd;
    gamai = gama;
    
else

    disp(['"change resolution" by computing 2D SPLINE INTERPOLATION res=',num2str(res)])
    
    % size for the interpolated grid
    Ni = res*(N-1)+1; % new size in y
    Mi = res*(M-1)+1; % new size in x
=======
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
>>>>>>> ameda_v2

    % elemental spacing for the interpolated grid
    dx = diff(x(1,1:2))/res;
    dy = diff(y(1:2,1))/res;

    % interpolated grid
<<<<<<< HEAD
    [xi,yi] = meshgrid([0:Mi-1]*dx+min(x(:)),[0:Ni-1]*dy+min(y(:)));
=======
    [xi,yi] = meshgrid(min(x(:)):dx:(Mi-1)*dx+min(x(:)),...
            min(y(:)):dy:(Ni-1)*dy+min(y(:)));
>>>>>>> ameda_v2

    % Increase resolution of the mask
    maski = interp2(x,y,mask,xi,yi);
    maski(isnan(maski) | maski < 1) = 0;
    
    % initialize interp fields
    ui = zeros([Ni Mi P]);
    vi = ui;

    % Increase resolution of fields (interp2 with regular grid)
    for i=1:P
<<<<<<< HEAD
        disp([' step ',num2str(i)])
        ui(:,:,i) = interp2(x,y,squeeze(u(:,:,i)),xi,yi,'*spline');
        vi(:,:,i) = interp2(x,y,squeeze(v(:,:,i)),xi,yi,'*spline');
    end
    
    % interpolated 2d parameters
    DXi = Dx/res;
    bi  = xi*0 + b0*res;
    bxi = xi*0 + b0*res;
    Rdi   = xi*0 + Rd(1,1);
    gamai = Rdi ./ (Dx*deg)/res;
    
end

% Save non interpolated fields
save([path_out,'fields',runname],'x','y','mask','u','v','b','bx','Dx','Rd','gama','-v7.3')

% Save interpolated fields
x=xi;
y=yi;
mask=maski;
u=ui;
v=vi;
b=bi;
bx=bxi;
Dx=Dxi;
Rd=Rdi;
gama=gamai;
save([path_out,'fields_inter',runname],'x','y','mask','u','v','b','bx','Dx','Rd','gama','-v7.3')
=======
        ui(:,:,i) = interp2(x,y,squeeze(u(:,:,i)),xi,yi,'spline');
        vi(:,:,i) = interp2(x,y,squeeze(v(:,:,i)),xi,yi,'spline');
    end
    
end

disp(' ')
>>>>>>> ameda_v2

