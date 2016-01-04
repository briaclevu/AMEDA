function psi=compute_psi(lon,lat,mask,u,v,ci,cj)
% compute_psi.m
%  compute_psi(lon,lat,mask,u,v,ci,cj) computes the streamfunction (PSI) field
%  by spatially integrating the u- and v-component of velocity within the 
%  area around the detected eddy center (ci,cj) defined in eddy_dim.m.
%  The underlaying assumption is that in the presence of an eddy the
%  the velocity field is characterized by weak divergence, so that contours
%  of PSI are tangential to the velocity vectors.
%  To reduce the error associated with this assumption PSI is computed
%  integrating u and v along two different paths, and the two fields are
%  then averaged.
%
%  - u and v are NxM matrices of the two component of velocity;
%  - lon and lat converted to km_di and km_dj which are NxM matrices
%      of the longitudinal and latitudinal 
%      spacing in km between grid points. (They can vary along i and j)
%  - ci and cj are indices of the center in the lon,lat grid
%
%  OUTPUT:
%  - psi is the NxM matrix of the streamfunction field
%
%  See also Appendix in Nencioli et al. (2009) for further details.
%
%-------------------------
%   Ver Dec 2015 by B. LE VU
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

% build mask
mask(mask==0)=nan;
[N0,M0] = size(mask);
% Nil if nan
u(isnan(u))=0;
v(isnan(v))=0;

%% prepare distance matrice for the 4 domains -----------------------

%----------------------------------------
% convert lon and lat into km distance matrices
km_di=nan(size(lon,1),size(lon,2)-1);
km_dj=nan(size(lat,1)-1,size(lat,2));
for i=1:size(lon,1)
    km_di(i,:)=sw_dist(lon(i,:),lat(i,:),'km');
end
for i=1:size(lat,2)
    km_dj(:,i)=sw_dist(lon(:,i),lat(:,i),'km');
end
    
%----------------------------------------
% length for the 4 domains (NE,SE,NO,SO)
lx1=size(u(cj:end,:),1);
lx2=size(u(1:cj,:),1);
ly1=size(u(:,ci:end),2);
ly2=size(u(:,1:ci),2);

%----------------------------------------
% adjust km_di and km_dj for trapz integration starting at ci, cj
% (add first row/column of zeros)
di=[km_di(:,1:ci-1) zeros(lx1+lx2-1,1) km_di(:,ci:end)];
dj=[km_dj(1:cj-1,:); zeros(1,ly1+ly2-1); km_dj(cj:end,:)];

%% compute streamfunction -----------------------

%----------------------------------------
% integrate first row of v along longitude (first term of eq.A2)
cx1=cumtrapz(v(cj,ci:end)).*di(cj,ci:end); % trapezoidal sum
cx2=-cumtrapz(v(cj,ci:-1:1)).*di(cj,ci:-1:1); % trapezoidal sum

% integrate first column of u along latitude (second term of eq.A3)
cy1=-cumtrapz(u(cj:end,ci)).*dj(cj:end,ci);
cy2=cumtrapz(u(cj:-1:1,ci)).*dj(cj:-1:1,ci);

%----------------------------------------
% expand the vectors into matrices to compute PSI
mcx11=cx1(ones(lx1,1),:);
mcx12=cx1(ones(lx2,1),:);
mcx21=cx2(ones(lx1,1),:);
mcx22=cx2(ones(lx2,1),:);
mcy11=cy1(:,ones(1,ly1));
mcy12=cy1(:,ones(1,ly2));
mcy21=cy2(:,ones(1,ly1));
mcy22=cy2(:,ones(1,ly2));

%----------------------------------------
% PSI from integrating v firts and then u
psi_xy11=(mcx11 - cumtrapz(u(cj:end,ci:end)).*dj(cj:end,ci:end)); %(eq. A2)
psi_xy12=(mcx12 + cumtrapz(u(cj:-1:1,ci:end)).*dj(cj:-1:1,ci:end)); %(eq. A2)
psi_xy21=(mcx21 - cumtrapz(u(cj:end,ci:-1:1)).*dj(cj:end,ci:-1:1)); %(eq. A2)
psi_xy22=(mcx22 + cumtrapz(u(cj:-1:1,ci:-1:1)).*dj(cj:-1:1,ci:-1:1)); %(eq. A2)

% concatene 4 parts (NE,SE,NO,SO)
psi_xy=[psi_xy22(end:-1:2,end:-1:2) psi_xy12(end:-1:2,:);...
	psi_xy21(:,end:-1:2) psi_xy11];

%----------------------------------------
% PSI from integrating u first and then v 
psi_yx11=(mcy11 + cumtrapz(v(cj:end,ci:end),2).*di(cj:end,ci:end)); %(eq. A3)
psi_yx21=(mcy12 - cumtrapz(v(cj:end,ci:-1:1),2).*di(cj:end,ci:-1:1)); %(eq. A3)
psi_yx12=(mcy21 + cumtrapz(v(cj:-1:1,ci:end),2).*di(cj:-1:1,ci:end)); %(eq. A3)
psi_yx22=(mcy22 - cumtrapz(v(cj:-1:1,ci:-1:1),2).*di(cj:-1:1,ci:-1:1)); %(eq. A3)

% concatene 4 parts (NE,SE,NO,SO)
psi_yx=[psi_yx22(end:-1:2,end:-1:2) psi_yx12(end:-1:2,:);...
        psi_yx21(:,end:-1:2) psi_yx11];

%----------------------------------------
% computed PSI as average between the two
psi=(psi_xy+psi_yx)/2.*mask;

%% Enlarge PSI into land by 1 pixel -----------------------
% PSI in land is an average of its 9 neighbours

mask(isnan(mask))=0;

for i=1:N0
    for j=1:M0
        if mask(i,j)==0 && ...
                sum(sum(mask(max(i-1,1):min(i+1,N0),max(j-1,1):min(j+1,M0))))~=0
            
            psi1=psi(max(i-1,1):min(i+1,N0),max(j-1,1):min(j+1,M0));
            psi(i,j)=nanmean(psi1(:));
            
        end
    end
end

