function [xbary,ybary,z,a,b,alpha,lim]=compute_ellip(xy,grid_ll)
%[xbary,ybary,z,a,b,alpha,lim]=compute_ellip(xy {,grid_ll})
%
% Compute the barycentre of a closed polygon 2xN array xy where
%   x or lon is in xy(1,:)
%   y or lat is in xy(2,:)
%   grid_ll =1  (default) if the coordinates are in (longitude,latitude)
%           =0 if coordinates are in 'km'
%
% and fit an ellipse using fitellipse
%
% (xbary,ybary) barycenter from coordinates of the polygon
%
%-------------------------
%  Jan 2016 by B. LE VU
%-------------------------
%
%=========================

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

xbary = somme_x/lim;
ybary = somme_y/lim;

%----------------------------------------
% Ellipse fitting
if lim >= 5
    try
        %----------------------------------------
        % Coord = coordinates of the polygon
        if grid_ll

            xs(1) = xbary;
            ys(1) = ybary;

            % initialise
            coord = zeros(2,lim);

            for pt=1:lim

                xs(2) = xy(1,pt);
                ys(2) = xy(2,pt);

                %----------------------------------------
                % distances in km of every point from the barycenter
                coord(1,pt) = sign(diff(xs)) * sw_dist([ybary ybary],xs,'km');
                coord(2,pt) = sign(diff(ys)) * sw_dist(ys,[xbary xbary],'km');
            end

            [~, a, b, alpha] = fitellipse(coord,'linear');
            z = [xbary,ybary];
            
        else

            coord = [xy(1,1:lim);xy(2,1:lim)];
            [z, a, b, alpha] = fitellipse(coord,'linear');
            
        end

        if a<0 || b<0 || z(1)==0
            z = [NaN NaN];
            a = NaN;
            b = NaN;
            alpha = NaN;
        end
        
    catch erre
        z = [NaN NaN];
        a = NaN;
        b = NaN;
        alpha = NaN;
    end
    
else
    
    z = [NaN NaN];
    a = NaN;
    b = NaN;
    alpha = NaN;
    
end

