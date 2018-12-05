function [CD,xy,allines,rmax,velmax,tau,deta,nrho,large,warn,calcul] =...
    eddy_dim(u,v,ssh,mask,x,y,centers,ii,farea,Rdarea,bxarea,plo)
%[CD,xy,allines,rmax,velmax,tau,deta,nrho,large,warn,calcul] =...
%               eddy_dim(u,v,ssh,mask,x,y,centers,ii {,f,Rd,bx,plo})
%
%  Computes the shape of the eddy defined by the iith center
%
%  - u and v are the 2D u and v velocity field of the step
%  - ssh is the 2D ssh field of the step (can be = [])
%  - centers is the potential eddy center structure
%  - x and y are the grid coordinates in the all domain of any step
%  - ii is the indice of the main eddy center
%  - farea is the Coriolis parameter at the center
%  - Rdarea is the deformation radius where the shape is computed
%  - bxarea is used to define the area where the shape is computed
%  - plo is a debug mod to check result of max_curve on a plot.
%       plot==0 by default
%
%  OUTPUT:
%  - CD is the (x,y) centers of the eddy and double eddy
%  - xy is the array containing x (first row) and
%    y (second row) cordinates of the vertices defining the eddy shape
%  - allines contains the profil V-R, tau, eta and nrho
%  - rmax is the radius of the contour xy
%  - velmax is the maximum mean velocity along the contour xy
%  - tau is the minimum turnover time inside the eddy contour
%  - deta is the deformation of the contour in absolute value 
%  - nrho is the part of the contour with negative curvature
%  - large and warn are flags that give some informations on the
%    procedure used to compute the eddy shape. See 'mod_eddy_shapes.m' for
%    further details.

%
% 'compute_psi' is used to compute the streamfunction field integrating u-
% and v-component in the geostrophic dynamic. Check the documentation in
% 'compute_psi.m' for further details.
%
% 'max_curve' is used to compute eddy shape (defined as the largest closed 
% contour of PSI around the eddy center across which velocity magnitude 
% increases). Check the documentation in 'max_curve.m' for further details.
%
%-------------------------
%   Ver. 3.2 Apr.2015 Briac Le Vu
%   Ver. 3.1 2014 LMD from Nencioli et al. routines
%-------------------------
%
%=========================

% load key_source and parameters
%----------------------------------------------
load('param_eddy_tracking')

% replace parameters by arguments
%----------------------------------------
bx = bxarea;
Rd = Rdarea;
f  = abs(farea);

if nargin<=11
    plo = 0;
end

%-----------------------------------------------------------
% main center type
type_c = centers.type(ii);

% main center coordinates
xy_cj = centers.y(ii);
xy_ci = centers.x(ii);

% indice of the main center (C_J and C_I) 
C_I = centers.i(ii);
C_J = centers.j(ii);

% all centers coordinates and type
type = centers.type;
centers_y = centers.y;
centers_x = centers.x;

% resize coordinate and velocity matrix 
% (making sure not to go outside the domain)
y = y(max(C_J-bx,1):min(C_J+bx,size(y,1)), ...
        max(C_I-bx,1):min(C_I+bx,size(y,2)));
x = x(max(C_J-bx,1):min(C_J+bx,size(x,1)), ...
        max(C_I-bx,1):min(C_I+bx,size(x,2)));
mask = mask(max(C_J-bx,1):min(C_J+bx,size(mask,1)), ...
        max(C_I-bx,1):min(C_I+bx,size(mask,2)));
v = v(max(C_J-bx,1):min(C_J+bx,size(v,1)), ...
        max(C_I-bx,1):min(C_I+bx,size(v,2)));
u = u(max(C_J-bx,1):min(C_J+bx,size(u,1)), ...
        max(C_I-bx,1):min(C_I+bx,size(u,2)));
if type_detection==2 || type_detection==3
    ssh = ssh(max(C_J-bx,1):min(C_J+bx,size(ssh,1)), ...
            max(C_I-bx,1):min(C_I+bx,size(ssh,2)));
    % determine contours to be scan
    Hs = floor(nanmin(ssh(:))):DH:ceil(nanmax(ssh(:)));
    % impose a limit to limit cpu time
    if length(Hs)>nH_lim
        Hs = nH_lim;
    end
end

% indice of the center in the small domain
[cj,ci] = find(y==xy_cj & x==xy_ci);

% coordinates of all eddy centers in the smaller area
% (contains at least xy_cj and xy_ci)
in = inpolygon(centers_x,centers_y,...
    [x(1,1) x(1,end) x(end,end) x(end,1)],...
    [y(1,1) y(1,end) y(end,end) y(end,1)]);
type_cts = type(in);
xy_ctsj = centers_y(in);
xy_ctsi = centers_x(in);

if type_detection==1 || type_detection==3
    
    % to calculate psi extrapole u and v to 0 in the land
    u(isnan(u)) = 0;
    v(isnan(v)) = 0;

    % compute streamfunction field from u,v in m/s/100 -> units ssh (m?)
    psi = compute_psi(x,y,mask,u*f/g*1e3,v*f/g*1e3,ci,cj,grid_ll);
    
    % determine contours to be scan
    H = floor(nanmin(psi(:))):DH:ceil(nanmax(psi(:)));
    
    % impose a limit to limit cpu time
    if length(H)>nH_lim
        H = nH_lim;
    end
    
end

%------------------------------------
% compute eddy shapes
switch type_detection
    
    case 1
        % compute eddy shape with psi
        [cd,eddy_lim,lines,rmax,velmax,tau,eta,nrho,large] =...
            max_curve(x,y,psi,xy_ci,xy_cj,type_cts,xy_ctsi,xy_ctsj,u,v,Rd,...
            H,n_min,k_vel_decay,nR_lim,Np,nrho_lim,grid_ll);
        calcul=1;
    
    case 2
        % compute eddy shape with ssh
        [cd,eddy_lim,lines,rmax,velmax,tau,eta,nrho,large] =...
            max_curve(x,y,ssh,xy_ci,xy_cj,type_cts,xy_ctsi,xy_ctsj,u,v,Rd,...
            Hs,n_min,k_vel_decay,nR_lim,Np,nrho_lim,grid_ll);
        calcul=0;

    case 3
        % compute eddy shape with ssh
        [cd,eddy_lim,lines,rmax,velmax,tau,eta,nrho,large] =...
            max_curve(x,y,ssh,xy_ci,xy_cj,type_cts,xy_ctsi,xy_ctsj,u,v,Rd,...
            Hs,n_min,k_vel_decay,nR_lim,Np,nrho_lim,grid_ll);
        calcul=0;

        if isnan(large(1))
            % compute eddy shape with psi
            [cd,eddy_lim,lines,rmax,velmax,tau,eta,nrho,large] =...
                max_curve(x,y,psi,xy_ci,xy_cj,type_cts,xy_ctsi,xy_ctsj,u,v,Rd,...
                H,n_min,k_vel_decay,nR_lim,Np,nrho_lim,grid_ll);
            calcul=1;
        end
end

%------------------------------------
% compute eddy features and flags
if ~calcul
    psi = ssh;
    H = Hs;
end

% deta initialisation
deta = nan(1,3);

% In case there is no contour with maximal velocity around the main center
% the eddy contour is associated to the last contour if any or erased
if isnan(large(1)) || isempty(eddy_lim{1})
    
    if isempty(eddy_lim{3})
        warn = 1;
        CD = [];
        xy = cell(1,3);
        allines = [];
    elseif rmax(3)<nR_lim*Rd && nrho(3)<nrho_lim
        eddy_lim{1} = eddy_lim{3};
        large(1) = 1;
        rmax(1) = rmax(3);
        velmax(1) = velmax(3);
        tau(1) = tau(3);
        eta(1) = eta(3);
        nrho(1) = nrho(3);
    end
    
end
    
if ~isempty(eddy_lim{3})

    warn = 0;
    CD = cd;
    xy = eddy_lim;
    allines = lines;
    
    % compute deta3
    in_eddy = inpolygon(x,y,xy{3}(1,:),xy{3}(2,:));
    if type_c == -1 % anticylone case
        deta(3) = max(psi(in_eddy~=0))-eta(3);
    elseif type_c == 1 % cyclone case
        deta(3) = min(psi(in_eddy~=0))-eta(3);
    end

    % compute deta1
    if ~isempty(eddy_lim{1})
        in_eddy = inpolygon(x,y,xy{1}(1,:),xy{1}(2,:));
        if type_c == -1 % anticylone case
            deta(1) = max(psi(in_eddy~=0))-eta(1);
        elseif type_c == 1 % cyclone case
            deta(1) = min(psi(in_eddy~=0))-eta(1);
        end
    end
    
    % compute deta2
    if ~isempty(eddy_lim{2})
        in_eddy = inpolygon(x,y,xy{2}(1,:),xy{2}(2,:));
        if type_c == -1 % anticylone case
            deta(2) = max(psi(in_eddy~=0))-eta(2);
        elseif type_c == 1 % cyclone case
            deta(2) = min(psi(in_eddy~=0))-eta(2);
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the shapes of the eddies detected in the domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plo && (~isempty(eddy_lim{1}) || ~isempty(eddy_lim{2}) || ~isempty(eddy_lim{3}))
    
 figure
 contour(x,y,psi,H)
 
 hold on
 quiver(x,y,u,v,'k')
 
 plot(xy_ctsi,xy_ctsj,'r*')
 plot(xy_ci,xy_cj,'k*')
 
 if ~isempty(xy{1})
    if large(1)==1
        plot(xy{1}(1,:),xy{1}(2,:),'.k');
    else
        plot(xy{1}(1,:),xy{1}(2,:),'-k','linewidth',3);
    end
 end
 
 if ~isempty(xy{2})
    plot(CD(1,:),CD(2,:),'bo')
    if large(2)==1
        plot(xy{2}(1,:),xy{2}(2,:),'.g');
    else
        plot(xy{2}(1,:),xy{2}(2,:),'-g','linewidth',3);
    end
 end
 
 if ~isempty(xy{3})
    plot(xy{3}(1,:),xy{3}(2,:),':b','linewidth',3);
 end
 
 hold off
 title(['Eddy ',num2str(ii)])
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
