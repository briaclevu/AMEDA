function psi = compute_psi(x,y,mask,u,v,ci,cj,grid_ll)
%[psi] = compute_psi(x,y,mask,u,v,ci,cj {,grid_ll} )
%
%  Computes the streamfunction (PSI) field by spatially integrating of
%  the u- and v-component of velocity in the geostrophic equilibrium
%  within the 4 area around a detected eddy center (ci,cj).
%  The underlaying assumption is that in the presence of an eddy the
%  velocity field is characterized by weak divergence, so that contours
%  of PSI are tangential to the velocity vectors.
%  To reduce the error associated with this assumption PSI is computed
%  integrating u and v along two different paths, and the two fields are
%  then averaged.
%
%  - x and y are NxM matrices of the coordinates
%  - u and v are NxM matrices of the two component of velocity;
%  - ci and cj are indices of the center in the grid
%  - grid_ll =1  (default) if the coordinates are in (longitude,latitude)
%            =0 if coordinates are in 'km'
%
%  OUTPUT:
%  - psi is the NxM matrix of the streamfunction field
%
%-------------------------
%   Jan 2016 by B. LE VU
%-------------------------
%
%=========================

% Default grid is (lon,lat)
if nargin==7
    grid_ll = 1;
end

% build mask
mask(mask==0) = nan;
[N0,M0] = size(mask);

% Nil if nan
u(isnan(u)) = 0;
v(isnan(v)) = 0;

%% prepare distance matrice for the 4 domains -----------------------

%----------------------------------------
% convert lon and lat into km distance matrices (can vary along i and j)
km_di = nan(size(x,1),size(x,2)-1);
km_dj = nan(size(y,1)-1,size(y,2));

for i=1:size(x,1)
    if grid_ll
        km_di(i,:) = sw_dist2(y(i,:),x(i,:));
    else
        km_di(i,:) = sqrt(diff(x(i,:)).^2 + diff(y(i,:)).^2); % km
    end
end

for i=1:size(y,2)
    if grid_ll
        km_dj(:,i) = sw_dist2(y(:,i),x(:,i));
    else
        km_dj(:,i) = sqrt(diff(x(:,i)).^2 + diff(y(:,i)).^2); % km
    end
end
    
%----------------------------------------
% length for the 4 domains (NE,SE,NO,SO)
lx1 = size(u(cj:end,:),1);
lx2 = size(u(1:cj,:),1);
ly1 = size(u(:,ci:end),2);
ly2 = size(u(:,1:ci),2);

%----------------------------------------
% adjust km_di and km_dj for trapz integration starting at ci, cj
% (add first row/column of zeros)
di = [km_di(:,1:ci-1) zeros(lx1+lx2-1,1) km_di(:,ci:end)];
dj = [km_dj(1:cj-1,:); zeros(1,ly1+ly2-1); km_dj(cj:end,:)];

%% compute streamfunction -----------------------

%----------------------------------------
% integrate first row of v along x
% (first term of the 4 parts of equations A1
cx1 =  cumtrapz(v(cj,ci:end))  .* di(cj,ci:end); % trapezoidal sum
cx2 = -cumtrapz(v(cj,ci:-1:1)) .* di(cj,ci:-1:1); % trapezoidal sum

% integrate first column of u along y
% (first term of the 4 parts of equations A2
cy1 = -cumtrapz(u(cj:end,ci))  .* dj(cj:end,ci);
cy2 =  cumtrapz(u(cj:-1:1,ci)) .* dj(cj:-1:1,ci);

%----------------------------------------
% expand the vectors into matrices to compute PSI
mcx11 = cx1(ones(lx1,1),:);
mcx12 = cx1(ones(lx2,1),:);
mcx21 = cx2(ones(lx1,1),:);
mcx22 = cx2(ones(lx2,1),:);
mcy11 = cy1(:,ones(1,ly1));
mcy12 = cy1(:,ones(1,ly2));
mcy21 = cy2(:,ones(1,ly1));
mcy22 = cy2(:,ones(1,ly2));

%----------------------------------------
% PSI from integrating v first and then u (4 parts of eq. A1)
psi_xy11 = (mcx11 - cumtrapz(u(cj:end,ci:end))   .* dj(cj:end,ci:end));
psi_xy12 = (mcx12 + cumtrapz(u(cj:-1:1,ci:end))  .* dj(cj:-1:1,ci:end));
psi_xy21 = (mcx21 - cumtrapz(u(cj:end,ci:-1:1))  .* dj(cj:end,ci:-1:1));
psi_xy22 = (mcx22 + cumtrapz(u(cj:-1:1,ci:-1:1)) .* dj(cj:-1:1,ci:-1:1));

% concatene 4 parts (NE,SE,NO,SO) (eq. A1)
psi_xy = [psi_xy22(end:-1:2,end:-1:2) psi_xy12(end:-1:2,:);...
        psi_xy21(:,end:-1:2) psi_xy11];

%----------------------------------------
% PSI from integrating u first and then v (4 parts of eq. A2)
psi_yx11 = (mcy11 + cumtrapz(v(cj:end,ci:end),2)   .* di(cj:end,ci:end));
psi_yx21 = (mcy12 - cumtrapz(v(cj:end,ci:-1:1),2)  .* di(cj:end,ci:-1:1));
psi_yx12 = (mcy21 + cumtrapz(v(cj:-1:1,ci:end),2)  .* di(cj:-1:1,ci:end));
psi_yx22 = (mcy22 - cumtrapz(v(cj:-1:1,ci:-1:1),2) .* di(cj:-1:1,ci:-1:1));

% concatene 4 parts (NE,SE,NO,SO) (eq. A2)
psi_yx = [psi_yx22(end:-1:2,end:-1:2) psi_yx12(end:-1:2,:);...
        psi_yx21(:,end:-1:2) psi_yx11];

%----------------------------------------
% computed PSI as average between the two (eq. A = (A1 + A2) /2)
psi = (psi_xy+psi_yx)/2.*mask;

%% Enlarge PSI into land by 1 pixel -----------------------

% PSI in first coastal pixel is an average of its 9 neighbours

mask(isnan(mask)) = 0;

for i=1:N0
    for j=1:M0
        if mask(i,j)==0 && ...
                sum(sum(mask(max(i-1,1):min(i+1,N0),max(j-1,1):min(j+1,M0))))~=0
            
            psi1     = psi(max(i-1,1):min(i+1,N0),max(j-1,1):min(j+1,M0));
            psi(i,j) = nanmean(psi1(:));
        end
    end
end

