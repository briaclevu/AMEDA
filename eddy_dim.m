function [CD,xy,allines,velmax,tau,deta,nrho,large,warn,calcul] =...
    eddy_dim(u,v,ssh,mask,x,y,centers,ii,bxarea,plo)
%[CD,xy,allines,velmax,tau,deta,nrho,large,warn,calcul] =...
%               eddy_dim(u,v,ssh,mask,x,y,centers,ii {,bx,plo})
%
%  Computes the shape of the eddy defined by the iith center
%
%  - u and v are the 2D u and v velocity field of the step
%  - ssh is the 2D ssh field of the step (can be = [])
%  - centers is the potential eddy center structure
%  - x and y are the grid coordinates in the all domain of any step
%  - bxarea is used to define the area where the shape is computed
%  - ii is the indice of the main eddy center
%  - plo is a debug mod to check result of max_curve on a plot.
%       plot==0 by default
%
%  OUTPUT:
%  - CD is the (x,y) centers of the eddy and double eddy
%  - xy is the array containing x (first row) and
%    y (second row) cordinates of the vertices defining the eddy shape
%  - allines are the profil 
%  - velmax is the maximum mean velocity along the contour xy
%  - tau is the minimum turnover time inside the eddy contour
%  - deta is the deformation of the contour in absolute value 
%  - nrho is the part of the contour with negative curvature
%  - large and warn are flags that give some informations on the
%    procedure used to compute the eddy shape. See 'mod_eddy_shapes.m' for
%    further details.

%
% 'compute_psi' is used to compute the streamfunction field integrating u-
% and v-component f velocity. Check the documentation in 'compute_psi.m' 
% for further details.
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
if nargin==9
    plo = 0;
    bx = bxarea;
elseif nargin==10
    bx = bxarea;
end

% main center coordinates
%-----------------------------------------------------------
c_y = centers.y(ii);
c_x = centers.x(ii);

% all centers coordinates
centers_y = centers.y;
centers_x = centers.x;

% find the indice of the main center in the coarse domain (C_J and C_I) 
if resol==1
    C_I = centers.i(ii);
    C_J = centers.j(ii);
else
    dist = sum(bsxfun(@minus, cat(3,c_x,c_y), cat(3,x,y)).^2,3);
    [C_J,C_I] = find(dist==min(dist(:)),1);
end

% center coordinates in the coarse domain
xy_cj = y(C_J,C_I);
xy_ci = x(C_J,C_I);

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
end

% indice of the center in the small domain
[cj,ci] = find(y==xy_cj & x==xy_ci);

% coordinates of all eddy centers in the smaller area
% (contains at least xy_cj and xy_ci)
in = inpolygon(centers_x,centers_y,...
    [x(1,1) x(1,end) x(end,end) x(end,1)],...
    [y(1,1) y(1,end) y(end,end) y(end,1)]);
xy_ctsj = centers_y(in);
xy_ctsi = centers_x(in);

if type_detection==1 || type_detection==3
    % to calculate psi extrapole u and v to 0 in the land
    u(isnan(u)) = 0;
    v(isnan(v)) = 0;

    % compute streamfunction field from u,v in m/s/100 -> units ssh (m?)
    psi = compute_psi(x,y,mask,u/100,v/100,ci,cj,grid_ll);
end

% compute eddy shapes
%------------------------------------
switch type_detection
    
    case 1
        % compute eddy shape with psi
        [cd,eddy_lim,lines,velmax,tau,eta,nrho,large] =...
            max_curve(x,y,psi,c_x,c_y,xy_ctsi,xy_ctsj,u,v,Rd,...
            H,n_min,k_vel_decay,R_lim,nrho_lim,grid_ll);
        calcul=1;
    
    case 2
        % compute eddy shape with ssh
        [cd,eddy_lim,lines,velmax,tau,eta,nrho,large] =...
            max_curve(x,y,ssh,c_x,c_y,xy_ctsi,xy_ctsj,u,v,Rd,...
            H,n_min,k_vel_decay,R_lim,nrho_lim,grid_ll);
        calcul=0;

    case 3
        % compute eddy shape with ssh
        [cd,eddy_lim,lines,velmax,tau,eta,nrho,large] =...
            max_curve(x,y,ssh,c_x,c_y,xy_ctsi,xy_ctsj,u,v,Rd,...
            H,n_min,k_vel_decay,R_lim,nrho_lim,grid_ll);
        calcul=0;

        if isnan(large(1))
            % compute eddy shape with psi
            [cd,eddy_lim,lines,velmax,tau,eta,nrho,large] =...
                max_curve(x,y,psi,c_x,c_y,xy_ctsi,xy_ctsj,u,v,Rd,...
                H,n_min,k_vel_decay,R_lim,nrho_lim,grid_ll);
            calcul=1;
        end
end

% compute eddy features and flags
%------------------------------------
if ~calcul
    psi = ssh;
end

% Rmove streamfunction curve closed around an island
% still a probleme because don't recover a smaller curve without island!!!
%in_eddy=find(inpolygon(x,y,eddy_lim{1}(1,:),eddy_lim{1}(2,:)));
%if min(mask(in_eddy))==0
%    large(1) = NaN;
%end

% deta initialisation
deta = nan(1,3);

% In case there is no streamfunction curve closed around the main center
% the eddy contour is erased
if isnan(large(1)) || isempty(eddy_lim{1})
    
    warn = 1;
    CD = [];
    xy = cell(1,3);
    allines = [];
    
else
    warn = 0;
    CD = cd;
    xy = eddy_lim;
    allines = lines;
    
    % compute deta1
    in_eddy = inpolygon(x,y,xy{1}(1,:),xy{1}(2,:));
    if max(psi(in_eddy~=0)) > eta(1)
        deta(1) = max(psi(in_eddy~=0))-eta(1);
    elseif min(psi(in_eddy~=0)) < eta(1)
        deta(1) = min(psi(in_eddy~=0))-eta(1);
    end
    
    % compute deta3
    in_eddy = inpolygon(x,y,xy{3}(1,:),xy{3}(2,:));
    if max(psi(in_eddy~=0)) > eta(3)
        deta(3) = max(psi(in_eddy~=0))-eta(3);
    elseif min(psi(in_eddy~=0)) < eta(3)
        deta(3) = min(psi(in_eddy~=0))-eta(3);
    end
    
    if ~isnan(large(2))
        
        % compute deta2
        in_eddy = inpolygon(x,y,xy{2}(1,:),xy{2}(2,:));
        if max(psi(in_eddy~=0)) > eta(2)
            deta(2) = max(psi(in_eddy~=0))-eta(2);
        elseif min(psi(in_eddy~=0)) < eta(2)
            deta(2) = min(psi(in_eddy~=0))-eta(2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the shapes of the eddies detected in the domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plo && (~isempty(eddy_lim{1}) || ~isempty(eddy_lim{2}))
    
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
        plot(xy{1}(1,:),xy{1}(2,:),'-k','linewidth',2);
    end
    plot(xy{3}(1,:),xy{3}(2,:),':k','linewidth',2);
 end
 
 if ~isempty(xy{2})
    plot(CD(1,:),CD(2,:),'bo')
    if large(2)==1
        plot(xy{2}(1,:),xy{2}(2,:),'.k');
    else
        plot(xy{2}(1,:),xy{2}(2,:),'-k','linewidth',2);
    end
 end
 
 hold off
 title(['Eddy ',num2str(ii)])
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
