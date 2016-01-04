function [vel]=integrate_vel(lon,lat,u,v,xdata,ydata)
% integrate_vel.m
%  integrate_vel(d_i,d_j,u,v,xdata,ydata) interpolate velocity
%  on a curve and integrate the projection along the curve.
%  This integration is used to check the increase in
%  velocity along the closed contour and an be saved as Velmax
%
%  - lon and lat are coordinate matrices
%  - u and v are  is the velocity fields, .
%  - xdata and ydata are the location coordinate of the curve
%      vertices under consideration;
%
%  OUTPUT:
%  - vel contains the velocity magnitude on the curve vertice
%
%-------------------------
%   Ver. dec.2015 B. LE VU
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

%-------------------------------------------------------------------
% Compute velocities fields Ui and Vi on the vertices by interpolate
% u and v fields on the contour xdata,ydata

if min(lon(1,:)==lon(end,:)) && min(lat(:,1)==lat(:,end))
    
    % regular grid (mesh). Simple interp2 fonction
    
    Ui=interp2(lon,lat,u,xdata(1:end-1),ydata(1:end-1));
    Vi=interp2(lon,lat,v,xdata(1:end-1),ydata(1:end-1));
    
else
    
    % irregular grid (not a mesh). loop on grid is necessary
    Ui=zeros(1,length(xdata)-1);
    Vi=zeros(1,length(xdata)-1);
    
    for j=1:length(xdata)-1
        
        % Find i and j index of closest grid point to vertices
        aa=find(lon<=xdata(j) & lat<=ydata(j));
        minD=sqrt((lon(aa)-xdata(j)).^2+(lat(aa)-ydata(j)).^2);
        [J,I]=ind2sub(size(lon),aa(minD==min(minD)));

        % Closest point position in fraction of index
        if I==size(lon,2) % To fix case where contour is out of domain 
            delta_i=1;
            I=I-1;
        else
            delta_i=(xdata(j)-lon(J,I))/(lon(J,I+1)-lon(J,I));
        end
        if J==size(lat,1) % To fix case where contour is out of domain
            delta_j=1;
            J=J-1;
        else
            delta_j=(ydata(j)-lat(J,I))/(lat(J+1,I)-lat(J,I));
        end
        % U and V fields on vertices (bi-linear interpolator)
        Ui(j)=u(J,I)*(1-delta_i)*(1-delta_j)+ ...
            u(J,I+1)*delta_i*(1-delta_j)+ ...
            u(J+1,I)*(1-delta_i)*delta_j+ ...
            u(J+1,I+1)*delta_i*delta_j;
        Vi(j)=v(J,I)*(1-delta_i)*(1-delta_j)+ ...
            v(J,I+1)*delta_i*(1-delta_j)+ ...
            v(J+1,I)*(1-delta_i)*delta_j+ ...
            v(J+1,I+1)*delta_i*delta_j;
    end
end

%-------------------------------------------------------------------
% Integrate u and v fields along a contour line by projecting (Ui,Vi)
% on a tangente define by its theta angle.
% Theta is computed giving more weight to the shortest edge
% where the phase is closer to the angle (Ui,Vi)

% Distance (D) and angle phase (P) for each edge + repeat the first
[swd,phase] = sw_dist(ydata([end-1,1:end]),xdata([end-1,1:end]),'km');

% Correcting term to resolve the right clockwise successive angle
cor = fix(abs(diff(phase))/180)*360;

% Tangente angle (theta) weighted by edge length (D1*P2 + D2*P1) / (D1 + D2)
theta = ( swd(1:end-1).*phase(2:end) +...
        swd(2:end).* ( phase(1:end-1) + cor ) ) ./...
        ( swd(1:end-1) + swd(2:end) );

% Compute integrated velocity
vel = abs( sum( dot([Ui;Vi],[cosd(theta);sind(theta)]) .*...
        (swd(1:end-1) + swd(2:end)) )/2/ sum(swd(2:end)));

