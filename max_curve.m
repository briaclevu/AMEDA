function [cd,eddy_lim,lines,velmax,tau,eta,large]=...
            max_curve(lon,lat,psi,ll_ci,ll_cj,ll_ctsi,ll_ctsj,u,v)
%  max_curve(lon,lat,psi,ll_ci,ll_cj,ll_ctsi,ll_ctsj,vel) computes the eddy shape
%  defined as the largest closed contour of the streamfunction (PSI) field in meter,
%  across which velocity magnitude increases.
%  The increase of velocity is checked between successive fraction (mm by mm)
%  of the closed contours of the field including the center. 2 closed
%  contour around the center is a minimum to properly check the eddy.
%  Velocity is defined as the mean value of each point of the contour
%  defined by the contourc function.
%  
%  If velocity never increases betwwen successive closed contour,
%  the shape is defined by simply the largest closed contour around the center.
%
%  - lon and lat are coordinates in a regular or irregular grid
%  - psi is the PSI field, from where the eddy shape is derived
%  - ll_ci and ll_cj are the coordinates of the detected eddy cente
%	   found by routine mod_eddy_centers in the small area
%  - ll_ctsi and ll_ctsj are coordinates of all potential centers
%	   in the small area
%  - vel is the velocity magnitude field, used to check the increase in
%          velocity across the closed contour.
%
%  OUTPUT:
%  - cd(1:2) are the centers coordinate  when 2 centers are recorded
%  - eddy_lim(1:2) are a 2xn array containing the position
%         of the n vertices that define eddy shape with 1 and with 2 centers;
%         (first row x, second y positions)
%  - lines record features of every streamlines scanned (1:number of centers
%    included; 2:ssh level; 3:mean radius; 4:mean velocity; 5:turnover time) 
%  - velmax(1:2) are the maximum mean velocity between n vertices along
%         the eddy_lim1 and 2
%  - tau(1:2) are minimum turnover time inside the eddy contour
%  - eta(1:2) are the ssh value for the eddy_lim1 and 2
%  - large(1:2) are flags for eddy_lim1 and 2 for largest contour
%    (1 if the contour is the largest and no "true" maximum is found)
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

% H is the scale of the PSI field to be determined
%H=double(floor(nanmin(psi(:))*1000):ceil(nanmax(psi(:))*1000))/1000;
global H
global streamlines

%-----------------------------------------------------------
% compute contourlines of the streamfunction psi (every 1mm)
% C is a vector containing all the coordinates of PSI contours

% test the regularity of the grid
if min(lon(1,:)==lon(end,:)) && min(lat(:,1)==lat(:,end))
    C = contourc(lon(1,:),lat(:,1),double(psi),H);
else
    figure('visible','off')
    C = contour(lon,lat,psi,H);
end

%-----------------------------------------------------------
% rearrange all the contours in C to the structure array 'isolines'
% sort by maximum latitude. Each element of isolines contains
% all the vertices of a given contour level of PSI

% fill the structure 'isolines'
isolines = struct('x',{},'y',{},'l',{});

% begin two counters
i = 1;
ii = 1;

while i < size(C,2)
    npoints = C(2,i);
    lvl(ii) = C(1,i);
    isolines(ii).x = C(1,i+1:i+npoints); % vertex lon's
    isolines(ii).y = C(2,i+1:i+npoints); % vertex lat's
    isolines(ii).l = max(C(2,i+1:i+npoints)); % max lat of a curve
    ii=ii+1;
    i=i+npoints+1;
end

% sort the contours according to their maximum latitude; this way the first
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
eddy_lim = cell(1,2); % lim with 1 and 2 centers
lines = []; % all streamlines features
large = nan(1,2);
velmax = zeros(1,2);
tau = nan(1,2);
eta = nan(1,2);

% to start properly fix the mean velocity and the eddy turnover time
Vmax = 0;
Tmin = [9999 9999];

% start the counter of isolines
i = 1;

%-----------------------------------------------------------
% scan all isolines until the eddy shape is determined
% and test the averaged velocity along each streamline
while i<=length(isolines)

    xdata = isolines(i).x; % vertex lon's
    ydata = isolines(i).y; % vertex lat's

    % conditions to have determine a contour
    % (isolines already sorted by maximum latitude)
    % 1) closed contours
    % 2) detected the right eddy center inside the polygon
    % 3) at least 1 contour without any other center
    % 4) no more than 3 centers
    % 5) record one contour with the right eddy center
    % 6) record the "true" maximum or by default the largest

    if xdata(1)==xdata(end) && ydata(1)==ydata(end) && ...
            inpolygon(ll_ci,ll_cj,xdata,ydata)
     
        % searchs the coordinates of the centers in the contour 
        IN = inpolygon(ll_ctsi,ll_ctsj,xdata,ydata);
        [~,p] = find(IN==1); % index of the coordinates in the contour
        nc = length(p); % number of center in the contour

        % test if the contour contains less than 2 centers
        if nc<=2

            % projection of velocities fields on a contour
            % for ir/regular grid and integrate these fields
            % along this contour to get the averaged velocity
            V = integrate_vel(lon,lat,u,v,xdata,ydata);
            
            % compute the perimeter
            P = sum(sw_dist(ydata(2,:),xdata(1,:),'km'));
            
            % compute the revolution time (Tau)
            T = P*1000/V/3600/24; % in days if vel is m/s

            % record every streamlines features
            if streamlines==1
                % Compute the R circle of a similar area
                R = rayon_moyen([xdata;ydata]);
                lines = [lines;nc lvl(i) R(1) V T];
            else
                lines = [lines;1];
            end

            % in case of 2 centers record (lon,lat) of the centers
            if nc==2
                cd = [ll_ctsi(p);ll_ctsj(p)];
            end
            
            % no closed contour met yet
            if Vmax==0
                
                % test if the first contour contains only 1 center
                if nc==1
                    % fix the test values
                    Vmax = V; % first value of velmax
                    Tmin(1) = min(Tmin(1),T); % first value of Tau
                else
                    i = length(isolines); % stop the scan
                end
                
            % closed contour already met and velocity is increasing
            elseif V>Vmax

                % update test values
                Vmax = V; % replace the velmax
                Tmin(nc) = min(Tmin(nc),T); % replace the Tau
                % record others index
                linesmax = [xdata;ydata];
                etamax = lvl(i);

                % record bigger contour around the single or 2 centers
                if length(xdata) > 6
                    
                    % test if no "true" maximum met yet
                    if large(nc)~=0

                        % replace previous contour
                        large(nc) = 1; % largest contour
                        velmax(nc) = Vmax;
                        eddy_lim{nc} = linesmax; % save the shape
                        tau(nc) = Tmin(nc); % save the turnover time
                        eta(nc) = etamax; % save the ssh contour
                    end
                end

           % velocity is decreasing more than 5% !!! new parameter !!!
            elseif V<Vmax*0.95
                
                % test if the last contour is the largest
                if large(nc)==1
                    large(nc) = 0; % we found an eddy with nc centers
                    
                % test if Vmax is higher then the existing "true" maximum
                elseif Vmax>velmax(nc) && velmax(nc)~=0
                    
                        % replace previous contour
                        velmax(nc) = Vmax;
                        eddy_lim{nc} = linesmax; % save the shape
                        tau(nc) = Tmin(nc); % save the turnover time
                        eta(nc) = etamax; % save the ssh contour
                end
            end

        % the contour contains more than 2 centers
        else
            i = length(isolines); % stop the scan
        end
    end
    i = i+1; % increase the counter
end

