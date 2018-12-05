function [curve,err] = compute_best_fit(lines,rmax,velmax)
% [curve,err] = compute_best_fit(lines,rmax,velmax)
%
% compute the best fit over an eddy profil V/Vmax = f(R/Rmax)
% the curve to fit must be the form:
%
%   Vr = Ro.e( (1-Ro^alpha) / alpha )
%
% with  Vr = V / Vmax
%       Ro = R / Rmax
%       alpha is the eddy degree (alpha=2 for a gaussian) to be fitted
%
% the curve is maximal for Vr = 1 (V=Vmax) at Ro = 1 (R=Rmax)
%
% Output the alpha with error, R² and chi_square
%
%-------------------------
%    June 2016 Briac LV
%-------------------------
%
%=========================

% initialize
curve = [];
err = [];

% read profil lines computed from max_curve
nc = lines(:,1);
eta = lines(:,2);
rmoy = lines(:,3);
vel = lines(:,4);
tau = lines(:,5);

% take part of the profil with no gap between point
% skip the potential gap at the beginning
indx = find(diff(rmoy(2:end)) > 2*mean(diff(rmoy(2:end))),1)+1; 

if isempty(indx)
    indx=length(rmoy);
end

% need decreasing points
indx1 = find(rmoy(1:indx)<rmax); % point before rmax
indx2 = find(rmoy(1:indx)>rmax); % point after rmax

% choose profil with 1 center, long enough, with a decreasing part
if nc(indx)<2 && length(indx1)>6 && length(indx2)>0.25*length(indx1) &&...
    min(vel(indx2)) < 0.9*velmax && max(vel(indx2)) < velmax

    %----------------------------------------------------------
    % Curve fitting from streamlines scanning (stop the scan at 'indx')
    % of an eddy (rmax,velmax) with mean radius 'rmoy' and
    % integrated velocity 'vel' for each streamlines
    %----------------------------------------------------------

    % change variable
    x=rmoy(1:indx)/rmax; %[0-2]*rmax
    y=vel(1:indx)/velmax; %[0-1]*velmax

    if size(x,2)~=1
    	x=x';y=y';
    end

    % define model function
    g = fittype(@(a,x) x.*exp((1-x.^a)/a),'depen','y','indep','x','coeff','a');
    
    % compute a fitting
    [curve,err] = fit(x,y,g,'Startpoint',1);
    
end

