function [R,A,P,ll] = mean_radius(xy,grid_ll)
%[R,A,P,ll] = mean_radius(xy {,grid_ll} )
%
% Compute the mean radius (R), the perimeter (P), the surface area (A)
% of a closed contour (xy) with is barycenter (ll).
%
% 2 cases:
% - (x,y) cartesian coordinates in 'km' if grid_ll=0
% - (lon,lat) earth coordinates if grid_ll=1 (default)
%
%-------------------------
%  Sept 2015 by B. LE VU
%-------------------------
%
%=========================

%----------------------------------------
% Default grid is (lon,lat)
if nargin==1
    grid_ll = 1;
end

%----------------------------------------
% Size of the polygon
lim = size(xy,2)-1;

%----------------------------------------
% Barycenter computation
somme_x = sum(xy(1,1:lim));
somme_y = sum(xy(2,1:lim));

ll(1) = somme_x/lim;
ll(2) = somme_y/lim;

%----------------------------------------
% Initialisation
distance=zeros(1,lim+1);
aire=zeros(1,lim);
param=zeros(1,lim);

%----------------------------------------
% Distance = distance between barycentre and every point of the polygon
xs(1)=ll(1);
ys(1)=ll(2);

for point=1:lim+1
    xs(2)=xy(1,point);
    ys(2)=xy(2,point);

    if grid_ll
        distance(point) = sw_dist2(ys,xs,'km');
    else
        distance(point) = sqrt(diff(xs).^2+diff(ys).^2); % km
    end
end

%----------------------------------------
% Distance2 = distance between 2 consecutives points of the polygon
if grid_ll
    distance2 = sw_dist2(xy(2,:),xy(1,:),'km');
else
    distance2 = sqrt(diff(xy(1,:)).^2+diff(xy(2,:)).^2); % km
end

%----------------------------------------
% Perimeter of the polygon computation
P = sum(distance2(:));

%----------------------------------------
% Surface area of the polygon computation
for point=1:lim
    a=distance(point);
    b=distance(point+1);
    c=distance2(point);
    s=(a+b+c)/2;
    aire(point)=sqrt(s*(s-a)*(s-b)*(s-c));
end

A = sum(real(aire(:)));

%% Mean radius computation (4 different radius)

for point=1:lim
    rayonmoy=(distance(point)+distance(point+1))/2;
    param(point)=aire(point)*rayonmoy;
end

%----------------------------------------
% Radius of a circle with equivalent surface area as the polygon (best)
R(1) = sqrt(A/(4*atan(1)));

%----------------------------------------
% Radius weight by the surface area of each polygon portion
R(2) = sum(param(:))/real(A);

%----------------------------------------
% Max radius from barycenter to points of the polygon
R(3) = max(distance(1:lim));

%----------------------------------------
% Mean radius from barycenter to points of the polygon
R(4) = mean(distance(1:lim));
