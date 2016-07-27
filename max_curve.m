function [cd,eddy_lim,lines,velmax,tau,eta,nrho,large] =...
            max_curve(x,y,psi,xy_ci,xy_cj,xy_ctsi,xy_ctsj,u,v,Rd,grid_ll)
%[cd,eddy_lim,lines,velmax,tau,eta,nrho,large] =...
%           max_curve(x,y,psi,xy_ci,xy_cj,xy_ctsi,xy_ctsj,u,v,Rd {,grid_ll} )
%
%  Computes the eddy shape defined as the largest closed contour of the
%  streamfunction (PSI/SSH) field in meter, across which velocity magnitude
%  increases.
%  The increase of velocity is checked between successive fraction set by H
%  of the closed contours and including the center. 2 closed
%  contours around the center is a minimum to properly check the eddy.
%  Velocity is defined as the mean value of the projected velocity for
%  each point of the contour defined by the contourc function.

%  
%  If velocity never increases betwwen successive closed contour,
%  the shape is defined by simply the largest closed contour around the center.
%
%  - x and y are the coordinates, earth (lon,lat) either cartesian (x,y)
%       in a regular or irregular grid of the small domain
%  - psi is the PSI or SSH field in the small domain
%  - xy_ci and xy_cj are the coordinates of the main detected eddy center
%       found by routine mod_eddy_centers in the small area
%  - xy_ctsi and xy_ctsj are coordinates of all potential centers
%	    in the small area
%  - u and v are the 2D u and v velocity field in the small area used to
%       compute the circulation around the eddy edge with integrate_vel.m
%  - Rd is the deformation radius
%  - grid_ll =1  (default) if the coordinates are in (longitude,latitude)
%            =0 if coordinates are in 'km'
%
%  OUTPUT:
%  - cd are the 2x2 centers coordinates (only when 2 centers are recorded)
%    (column 1 for the main center and 2 for the second)
%    (row 1 for x and 2 for y)
%  - eddy_lim {1x3} are the 3 2xn array containing the position
%    of the n vertices that define the eddy shape with 1 and with 2 centers
%    (first row is x, second is y positions) and also the last contour with
%    1 center.
%  - lines are px5 recorded features of the p streamlines scanned. Features
%    following the row number:
%       1:number of centers included
%       2:ssh level
%       3:mean radius
%       4:mean velocity
%       5:turnover time
%  - velmax is the maximum  1x3 mean velocity between n vertices along
%    the eddy_lim1, eddy_lim2 and eddy_lim3
%  - tau is the 1x2 minimum turnover times inside the eddies contour
%  - eta is the 1x3 PSI/SSH values for the eddy_lim1 and eddy_lim2
%  - nrho is the 1x2 part of the contour with negative curvature
%  - large is the 1x2 flags for eddy_lim1 and 2 for largest contour
%    (1 if the contour is the largest and no "true" maximum is found)
%
%-------------------------
%   Ver. 3.2 Apr.2015 Briac Le Vu
%   Ver. 3.1 2014 LMD from nencioli et al. routines
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

global streamlines
global H

% Default grid is (lon,lat)
if nargin==10
    grid_ll = 1;
end

%-----------------------------------------------------------
% compute contourlines of the streamfunction psi (H contour)
% C is a vector containing all the coordinates of PSI contours

% test the regularity of the grid
if min(x(1,:)==x(end,:)) && min(y(:,1)==y(:,end))
    C = contourc(x(1,:),y(:,1),double(psi),H);
else
    figure('visible','off')
    C = contour(x,y,psi,H);
end

%-----------------------------------------------------------
% rearrange all the contours in C to the structure array 'isolines'
% sort by maximum y. Each element of isolines contains
% all the vertices of a given contour level of PSI

% fill the structure 'isolines'
isolines = struct('x',{},'y',{},'l',{});

% begin two counters
i = 1;
ii = 1;

while i < size(C,2)
    npoints = C(2,i);
    lvl(ii) = C(1,i);
    isolines(ii).x = C(1,i+1:i+npoints); % vertex x's
    isolines(ii).y = C(2,i+1:i+npoints); % vertex y's
    isolines(ii).l = max(C(2,i+1:i+npoints)); % max y of a curve
    ii=ii+1;
    i=i+npoints+1;
end

% sort the contours according to their maximum y; this way the first
% closed contour across which velocity increases will also be the largest
% one (it's the one which extend further north).
[~,order] = sort([isolines(:).l],'ascend');
isolines = isolines(order);
lvl = lvl(order);

%-----------------------------------------------------------
% intialize the variables
% fix eddy contour and property empty
% if no closed contour exists around the center
cd = nan(2); % centers in the closed contour
eddy_lim = cell(1,3); % lim with 1 and 2 centers
lines = []; % all streamlines features
large = nan(1,2);
velmax = zeros(1,3);
tau = nan;
eta = nan(1,3);
nrho = nan(1,2);

% to start properly fix the mean velocity and the eddy turnover time
Vmax = 0;
Tmin = 9999;

% start the counter of isolines
i = 1;

%-----------------------------------------------------------
% scan all isolines until the eddy shape is determined
% and test the averaged velocity along each streamline
while i<=length(isolines)

    xdata = isolines(i).x; % vertex x's
    ydata = isolines(i).y; % vertex y's

    % conditions to have determine a contour
    % (isolines already sorted by maximum y)
    % 1) closed contours
    % 2) detected the right eddy center inside the polygon
    % 3) at least 1 contour without any other center
    % 4) no more than 2 centers
    % 5) record two contours with the right eddy center ('max' and 'end')
    % 6) record the "true" maximum or by default the largest
    
    if xdata(1)==xdata(end) && ydata(1)==ydata(end) && ...
            inpolygon(xy_ci,xy_cj,xdata,ydata) && length(xdata) >= 4
     
        % searchs the coordinates of the centers in the contour 
        IN = inpolygon(xy_ctsi,xy_ctsj,xdata,ydata);
        [~,p] = find(IN==1); % index of the coordinates in the contour
        nc = length(p); % number of center in the contour

        % test if the contour contains less than 2 centers
        if nc<=2

            % projection of velocities fields on a contour
            % for ir/regular grid and integrate these fields
            % along this contour to get the averaged velocity
            V = integrate_vel(x,y,u,v,xdata,ydata,grid_ll);
            
            % compute the R circle of a similar area
            R = mean_radius([xdata;ydata],grid_ll);
            
            % compute the local curvature (C) and the segment length (P)
            [C,P] = compute_curve([xdata;ydata],Rd,grid_ll);
            
            % compute the revolution time (Tau)
            T = sum(P(1:end-1))*1000/V/3600/24; % in days if vel is m/s
            
            % part of the contour with negative curvature weighted by
            % the discrete curvature
            N = abs(sum(P(C<0).*C(C<0))/(2*pi));
            
            % record every streamlines features
            if streamlines
                lines = [lines;nc lvl(i) R(1) V T];
            else
                lines = [lines;1];
            end

            % in case of 1 center record the last shape
            if nc==1
                velmax(3) = V;
                eddy_lim{3} = [xdata;ydata]; % save the last shape
                eta(3) = lvl(i); % save the ssh contour
            % in case of 2 centers record (x,y) of the centers
            elseif nc==2
                cd = [xy_ctsi(p);xy_ctsj(p)];
            end
            
            % no closed contour met yet
            if Vmax==0
                
                % test if the first contour contains only 1 center
                % it is not too big or too much concave segment
                if nc==1 && R(1)<5*Rd && N<0.33
                    % fix the test values
                    Vmax = V; % first value of velmax
                    Tmin = min(Tmin,T); % first value of Tau
                    % record the first shape which can be the last
                    velmax(3) = V;
                    eddy_lim{3} = [xdata;ydata]; % save the last shape
                    eta(3) = lvl(i); % save the ssh contour
                else
                    return % stop the scan
                end
                
            % closed contour already met and velocity is increasing
            elseif V>Vmax

                % Indice of the contour
                I = i;
		        % update test values
                Vmax = V; % replace the velmax
                Tmin = min(Tmin,T); % replace the Tau
                % record others index
                linesmax = [xdata;ydata];
                etamax = lvl(i);
                rmax = R(1);
                nrhomax = N;

                % record bigger contour around the single or 2 centers
                % if no "true" maximum met yet
                if large(nc)~=0

                    % replace previous contour
                    large(nc) = 1; % largest contour
                    % record eddy{1} only if R<Rlim and N<Nlim
                    if nc==2 || rmax<5*Rd && nrhomax<0.33
                        velmax(nc) = Vmax;
                        eddy_lim{nc} = linesmax; % save the shape
                        tau = Tmin; % save the turnover time
                        eta(nc) = etamax; % save the ssh contour
                        nrho(nc) = nrhomax; % save the local curvature
                    end                            
                end

            % velocity is decreasing more than 5% !!! new parameter !!!
            elseif V<Vmax*0.95
                
                % test if the last contour is the largest and not the previous one
                if large(nc)==1 && i-I>1
                    
                    large(nc) = 0; % we found an eddy with nc centers
                    
                % test if Vmax is higher then the existing "true" maximum
                elseif Vmax>velmax(nc) && velmax(nc)~=0
                    % record eddy{1} only if R<Rlim and N<Nlim
                    if nc==2 || rmax<5*Rd && nrhomax<0.33
                        % replace previous contour
                        velmax(nc) = Vmax;
                        eddy_lim{nc} = linesmax; % save the shape
                        tau = Tmin; % save the turnover time
                        eta(nc) = etamax; % save the ssh contour
                        nrho(nc) = nrhomax; % save the local curvature
                    end
                end
            end

        % the contour contains more than 2 centers
        else
            return % stop the scan
        end
    end
    i = i+1; % increase the counter
end
