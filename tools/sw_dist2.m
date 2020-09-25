function [dist] = sw_dist2(lat,lon)

% SW_DIST    Distance between two lat,lon coordinates in km
%===================================================================
% SW_DIST  $Id: sw_dist.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
%          Copyright (C) CSIRO, Phil Morgan & Steve Rintoul 1992.
%
% USAGE:  [dist] = sw_dist(lat,lon)
%
% DESCRIPTION:
%   Calculate distance between two positions on glode using the "Plane
%   Sailing" method.  Also uses simple geometry to calculate the bearing of
%   the path between position pairs.
%
% INPUT:
%    lat      = decimal degrees (+ve N, -ve S) [- 90.. +90]
%    lon      = decimal degrees (+ve E, -ve W) [-180..+180]
%
% OUTPUT:
%    dist        = distance between positions in units
%
% AUTHOR:   Phil Morgan and Steve Rintoul 92-02-10
%           Simplified by B. LE VU (2016) for AMEDA purpose (2 points only)
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE:
%    The PLANE SAILING method as descriibed in "CELESTIAL NAVIGATION" 1989 by
%    Dr. P. Gormley. The Australian Antartic Division.
%==================================================================

% Modifications
% 99-06-25. Lindsay Pender, Function name change from distance to sw_dist.
% 99-06-25. Lindsay Pender, Fixed transpose of row vectors.

% CALLER:   general purpose
% CALLEE:   none

%-----------------
% DEFINE CONSTANTS
%-----------------
DEG2RAD = (2*pi/360);
DEG2NM  = 60;
NM2KM   = 1.8520;    % Defined in Pond & Pickard p303.

% BEGIN
npositions = length(lat);
ind=1:npositions-1;     % index to first of position pairs

dlon = diff(lon);
if any(abs(dlon)>180)
   flag = find(abs(dlon)>180);
   dlon(flag)= -sign(dlon(flag)) .* (360 - abs(dlon(flag)) );
end %if
latrad = abs(lat*DEG2RAD);
dep    = cos( (latrad(ind+1)+latrad(ind))./2 ) .* dlon;
dlat   = diff(lat);
dist   = DEG2NM*sqrt(dlat.^2 + dep.^2);  % in n.miles

dist = dist * NM2KM;% in km

return
