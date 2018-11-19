function [vel] = integrate_vel(x,y,u,v,xdata,ydata,grid_ll)
%[vel] = integrate_vel(x,y,u,v,xdata,ydata {,grid_ll} )
%
%  Interpolate velocity (u,v) on a polygon (xdata,ydata) and integrate
%  the projection along this polygon. The result of the integration
%  is used to check the increase in velocity for successive streamlines
%  in max_curve.m and is saved as 'vmax' and in 'lines' by max_curve.m
%
%  - x and y are coordinate matrices
%  - u and v are  is the velocity fields component
%  - xdata and ydata are the location coordinate of the polygon
%      vertices under consideration.
%  - grid_ll =1  (default) if the coordinates are in (longitude,latitude)
%            =0 if coordinates are in 'km'
%
%  OUTPUT:
%  - vel contains the velocity magnitude on the polygon vertices
%
%-------------------------
%   Ver. dec 2015 B. LE VU
%   Ver. 3.1 2014 LMD
%-------------------------
%
%=========================

% Default grid is (lon,lat)
if nargin==6
    grid_ll = 1;
end

%-------------------------------------------------------------------
% Compute velocities fields Ui and Vi on the vertices by interpolate
% u and v fields on the polygon vertices (xdata,ydata)

if min(x(1,:)==x(end,:)) && min(y(:,1)==y(:,end))
    
    % regular grid (mesh). Simple interp2 fonction
    Ui = interp2(x,y,u,xdata(1:end-1),ydata(1:end-1));
    Vi = interp2(x,y,v,xdata(1:end-1),ydata(1:end-1));
    
else
    
    % irregular grid (not a mesh). loop on grid is necessary
    Ui = zeros(1,length(xdata)-1);
    Vi = zeros(1,length(xdata)-1);
    
    for j=1:length(xdata)-1
        
        % Find i and j index of closest grid point to vertices
        aa = find(x<=xdata(j) & y<=ydata(j));
        minD = sqrt((x(aa)-xdata(j)).^2+(y(aa)-ydata(j)).^2);
        [J,I] = ind2sub(size(x),aa(minD==min(minD)));

        % Closest point position in fraction of index
        if I==size(x,2) % To fix case where contour is out of domain 
            delta_i = 1;
            I = I - 1;
        else
            delta_i = (xdata(j)-x(J,I)) / (x(J,I+1)-x(J,I));
        end
        if J==size(y,1) % To fix case where contour is out of domain
            delta_j = 1;
            J = J - 1;
        else
            delta_j = (ydata(j)-y(J,I)) / (y(J+1,I)-y(J,I));
        end
        % U and V fields on vertices (bi-linear interpolator)
        Ui(j) = u(J,I) * (1-delta_i) * (1-delta_j) + ...
            u(J,I+1) * delta_i * (1-delta_j) + ...
            u(J+1,I) * (1-delta_i) * delta_j + ...
            u(J+1,I+1) * delta_i * delta_j;
        Vi(j) = v(J,I) * (1-delta_i) * (1-delta_j) + ...
            v(J,I+1) * delta_i * (1-delta_j) + ...
            v(J+1,I) * (1-delta_i) * delta_j + ...
            v(J+1,I+1) * delta_i * delta_j;
    end
end

%-------------------------------------------------------------------
% Integrate u and v fields along a polygon line by projecting (Ui,Vi)
% on a tangente define by its theta angle.
% Theta is computed giving more weight to the shortest edge
% where the phase is closer to the angle (Ui,Vi)

% Distance (D) and angle phase (P) for each edge + repeat the first
if grid_ll
    [swd,phase] = sw_dist(ydata([end-1,1:end]),xdata([end-1,1:end]),'km');
else
    swd = sqrt(diff(ydata([end-1,1:end])).^2 +...
            diff(xdata([end-1,1:end])).^2); % km
    phase = angle(diff(xdata([end-1,1:end])) +...
            diff(ydata([end-1,1:end]))*sqrt(-1))...
            / (2*pi)*360; % deg
end

% Correcting term to resolve the right clockwise successive angle
cor = fix( abs(diff(phase))/180 ) * 360;

% Tangente angle (theta) weighted by edge length (D1*P2 + D2*P1) / (D1 + D2)
theta = ( swd(1:end-1) .* phase(2:end) +...
        swd(2:end) .* ( phase(1:end-1) + cor ) ) ./...
        ( swd(1:end-1) + swd(2:end) );

% Compute integrated velocity
vel = abs( sum( dot([Ui;Vi],[cosd(theta);sind(theta)]) .*...
        (swd(1:end-1) + swd(2:end)) )/2/ sum(swd(2:end)));

