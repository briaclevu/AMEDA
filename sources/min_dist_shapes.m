function ds=min_dist_shapes(xy1,xy2,grid_ll) 
%[ds] = min_dist_shapes(lonlat1,lonlat2 {,grid_ll} )
%
%  Calcul the minimal distance in 'km' between 2 eddies edges identified
%  by xy1 & xy2 their contour polygon coordinates from eddy_dim.m in :
%  - earth (lon,lat) if grid_ll=1 (default)
%  - cartesian (x,y) in 'km' if grid_ll=0
%
%-------------------------
%  Sept 2015 by B. LE VU
%-------------------------
%
%=========================

% Default grid is (lon,lat)
if nargin==2
    grid_ll=1;
end

% Initialisation
ds=nan(length(xy1),length(xy2));

% Double loop (didn't find anything better... any idea!?)
for i=1:length(xy1)
    for j=1:length(xy2)
        ys=[xy1(2,i),xy2(2,j)];
        xs=[xy1(1,i),xy2(1,j)];
        if grid_ll
            ds(i,j) = sw_dist2(ys,xs);
        else
            ds(i,j) = sqrt(diff(xs).^2+diff(ys).^2); % km
        end
    end
end

% Minimal distance
ds=min(ds(:));
