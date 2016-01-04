function [CD,lonlat,allines,velmax,tau,deta,large,warn,box,calcul]=...
    eddy_dim(u,v,ssh,mask,lon,lat,centers,fac,rd,ii)
%  eddy_dim(u,v,mask,lon,lat,centers,fac,rd,ii)
%  computes the shape of the eddy defined by the iith center
%
%  - u and v are the 2D u and v velocity field of the step
%  - ssh is the 2D ssh field of the step
%  - centers is the eddy center matrix
%  - lon and lat are longitudes and latitudes of points in the domain;
%  - fac is the factor to enlarge the area where the shape is computed;
%  - rd is used to define the area where the shape is computed;
%  - ii is the indice of the main eddy center
%
%  OUTPUT:
%  - CD is the (lon,lat) centers of the eddy and double eddy
%  - lonlat is the array containing longitudes (first row) and
%    latitudes (second row) of the vertices defining the eddy shape;
%  - velmax is the maximum mean velocity along the contour lonlat
%  - tau is the minimum turnover time inside the eddy contour
%  - deta is the deformation of the contour in absolute value 
%  - box is the flag which will be used in mod_eddy_shapes, that indicates if
%	   no maximum velocity has been met (1) and permit to increase the small area
%  - large, warn and box are flags that give some informations on the
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
%   Ver. 3.1 2014 LMD
%   Ver. 2.1 Oct.2012
%   Ver. 2.0 Jan.2012
%   Ver. 1.3 Apr.2011
%   Ver. 1.2 May.2010
%   Ver. 1.1 Dec.2009
%   Authors: Francesco Nencioli, francesco.nencioli@univ-amu.fr
%            Charles Dong, cdong@atmos.ucla.edu
%-------------------------
%
% Copyright (C) 2009-2012 Francesco Nencioli and Charles Dong
%
% This file is part of the Vector-Geometry Eddy Detection Algorithm.
%
% The Vector-Geometry Eddy Detection Algorithm is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
% 
% The Vector-Geometry Eddy Detection Algorithm is distributed in the 
% hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
% the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
% PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Vector-Geometry Eddy Detection Algorithm.  
% If not, see <http://www.gnu.org/licenses/>.
%
%=========================

global type_detection

% debug mode
plo = 0;

% increase the dimensions of the area where eddy shape is computed
rd = rd*fac;

%-----------------------------------------------------------
% main center coordinates
c_lat = centers.lat(ii);
c_lon = centers.lon(ii);

% all centers coordinates
centers_lat = centers.lat;
centers_lon = centers.lon;

% find the indice of the main center in the coarse domain (C_J and C_I) 
[~,C_J] = min(abs(lat(:,1)-centers.lat(ii)));
[~,C_I] = min(abs(lon(C_J,:)-centers.lon(ii)));

% center coordinates in the coarse inter domain
ll_cj = lat(C_J,C_I);
ll_ci = lon(C_J,C_I);

% resize coordinate and velocity matrix 
% (making sure not to go outside the domain)
lat = lat(max(C_J-rd,1):min(C_J+rd,size(lat,1)), ...
        max(C_I-rd,1):min(C_I+rd,size(lat,2)));
lon = lon(max(C_J-rd,1):min(C_J+rd,size(lon,1)), ...
        max(C_I-rd,1):min(C_I+rd,size(lon,2)));
mask = mask(max(C_J-rd,1):min(C_J+rd,size(mask,1)), ...
        max(C_I-rd,1):min(C_I+rd,size(mask,2)));
v = v(max(C_J-rd,1):min(C_J+rd,size(v,1)), ...
        max(C_I-rd,1):min(C_I+rd,size(v,2)));
u = u(max(C_J-rd,1):min(C_J+rd,size(u,1)), ...
        max(C_I-rd,1):min(C_I+rd,size(u,2)));
if type_detection==2 || type_detection==3
    ssh=ssh(max(C_J-rd,1):min(C_J+rd,size(ssh,1)), ...
        max(C_I-rd,1):min(C_I+rd,size(ssh,2)));
end

% indice of the center in the small domain
[cj,ci] = find(lat==ll_cj & lon==ll_ci);

% coordinates of all eddy centers in the smaller area
% (contains at least ll_cj and ll_ci)
in = inpolygon(centers_lon,centers_lat,...
    [lon(1,1) lon(1,end) lon(end,end) lon(end,1)],...
    [lat(1,1) lat(1,end) lat(end,end) lat(end,1)]);
ll_ctsj = centers_lat(in);
ll_ctsi = centers_lon(in);

if type_detection==1 || type_detection==3
    %------------------------------------
    % compute streamfunction field from u,v in m/s/100 -> units ssh (m?)
    psi = compute_psi(lon,lat,mask,u/100,v/100,ci,cj);
end

% compute eddy shapes
%------------------------------------
if type_detection==1
    %------------------------------------
    % compute eddy shape with psi
    [cd,eddy_lim,lines,velmax,tau,eta,large] =...
        max_curve(lon,lat,psi,c_lon,c_lat,ll_ctsi,ll_ctsj,u,v);
    calcul=1;
elseif type_detection==2
    %------------------------------------
    % compute eddy shape with ssh
    [cd,eddy_lim,lines,velmax,tau,eta,large] =...
        max_curve(lon,lat,ssh,c_lon,c_lat,ll_ctsi,ll_ctsj,u,v);
    calcul=0;
elseif type_detection==3
    %------------------------------------
    % compute eddy shape with ssh
    [cd,eddy_lim,lines,velmax,tau,eta,large] =...
        max_curve(lon,lat,ssh,c_lon,c_lat,ll_ctsi,ll_ctsj,u,v);
    calcul=0;
    if isnan(large(1))
        %------------------------------------
        % compute eddy shape with psi
        [cd,eddy_lim,lines,velmax,tau,eta,large] =...
            max_curve_inter_3(lon,lat,psi,c_lon,c_lat,ll_ctsi,ll_ctsj,u,v);
        calcul=1;
    end
end

% compute eddy features and flags
%------------------------------------
if ~calcul
    psi = ssh;
end

% deta initialisation
deta = nan(1,2);

% in case there is no streamfunction curve closed around the right center
% the eddy is erased
if isnan(large(1))
    warn = 1;
    box = 0;
    CD = [];
    lonlat = cell(1,2);
    allines = [];
else
    warn = 0;
    CD = cd;
    lonlat = eddy_lim;
    allines = lines;
    
    % enlarge the box if no "true" maximum found
    box = large(1);
    % compute deta1
    in_eddy = inpolygon(lon,lat,lonlat{1}(1,:),lonlat{1}(2,:));
    if max(psi(in_eddy~=0)) > eta(1)
        deta(1) = max(psi(in_eddy~=0))-eta(1);
    elseif min(psi(in_eddy~=0)) < eta(1)
        deta(1) = min(psi(in_eddy~=0))-eta(1);
    end
    
    if ~isnan(large(2))
        % enlarge the box if no "true" maximum found
        box = large(2);
        % compute deta2
        in_eddy = inpolygon(lon,lat,lonlat{2}(1,:),lonlat{2}(2,:));
        if max(psi(in_eddy~=0)) > eta(2)
            deta(2) = max(psi(in_eddy~=0))-eta(2);
        elseif min(psi(in_eddy~=0)) < eta(2)
            deta(2) = min(psi(in_eddy~=0))-eta(2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncomment to plot the shapes of the eddies detected in the domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plo && (~isempty(eddy_lim{1}) || ~isempty(eddy_lim{2}))
 global H
 figure
 contour(lon,lat,psi,H)
 hold on
 quiver(lon,lat,u,v,'k')
 plot(ll_ctsi,ll_ctsj,'r*')
 plot(ll_ci,ll_cj,'k*')
 if ~isempty(lonlat{1})
    if large(1)==1
        plot(lonlat{1}(1,:),lonlat{1}(2,:),'.k');
    else
        plot(lonlat{1}(1,:),lonlat{1}(2,:),'-k','linewidth',2);
    end
 end
 if ~isempty(lonlat{2})
    plot(CD(1,:),CD(2,:),'bo')
    if large(2)==1
        plot(lonlat{2}(1,:),lonlat{2}(2,:),'.k');
    else
        plot(lonlat{2}(1,:),lonlat{2}(2,:),'-k','linewidth',2);
    end
 end
 hold off
 title(['Eddy ',num2str(ii)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
