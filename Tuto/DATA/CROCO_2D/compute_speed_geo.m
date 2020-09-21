function [ugeo,vgeo] = compute_speed_geo(n,x,y,ssh)
% function [ugeo,vgeo] = compute_speed_geo(n,x,y,ssh)
%
%  compute geostrophique field by a n-stencil centered finite difference 
%  with non-uniform grid spacing in x and y from a ssh field
%
%-------------------------
%  April 2020 Briac Le Vu nased on https://github.com/mluhar/dynamicblade 
%-------------------------
%
%=========================

%----------------------------------------
% Calculate coriolis parameter
T = 3600*24; % a day period in seconds
f = 4*pi/T*sind(y); % in s-1

% define constants
earth_radius = 6378.137; % km

% kilometer (km) per degree of latitude
R = earth_radius*pi/180; % 111.320m

% gravitation constant
g = 9.8; % m.s-2

%----------------------------------------
% Calculation of finite spatial element
[~,dx] = gradient(x);
[dy,~] = gradient(y);

% Calcul finite spatial element in km
dcy=cumsum(dy,2)*R;
dcx=cumsum(dx*R.*cosd(y),1);

% in meters
dcx = dcx*1000; % m
dcy = dcy*1000; % m

%----------------------------------------
% initialisation
ugeo = nan(size(ssh));
vgeo = nan(size(ssh));

%----------------------------------------
disp('Compute geostrophic speed')
for t=1:size(ssh,3)

  disp(['  step ',num2str(t),'/',num2str(size(ssh,3))])

  % u from n-stencil centered difference of dy
  for i=1:size(dcy,1)
  
    % mask land and island along y(i,:)
    idnan = isnan(ssh(i,:,t));

    if length(dcy(i,~idnan))>n

      % the n dy Weight at every y(i,j)
      Dynan = fdmatrix(dcy(i,~idnan),1,n);
  
      % number of remaining ssh(i,:)
      Nnan = length(ssh(i,~idnan,t));

      % prepare the ssh matrix of ssh along y(i,:)
      sshi = repmat(ssh(i,~idnan,t),[Nnan 1]);

      % u(i,j) is the sum of Dy weighted ssh along y(i,:)
      ugeo(i,~idnan,t) = -g./f(i,~idnan).*nansum(Dynan.*sshi,2)';

    end

  end

  % v from n-stencil centered difference of dx latitude by latitude j
  for j=1:size(dcx,2)
  
    % mask land and island along x(:,j)
    idnan = isnan(ssh(:,j,t));

    if length(dcx(~idnan,j))>n

      % the n dx Weight at every x(i,j)
      Dxnan = fdmatrix(dcx(~idnan,j),1,n);
  
      % number of remaining ssh(:,j)
      Nnan = length(ssh(~idnan,j,t));
  
      % prepare the ssh matrix of ssh along x(:,j)
      sshj = repmat(ssh(~idnan,j,t)',[Nnan 1]);
  
      % v(i,j) is the sum of Dx weighted ssh along x(:,j)
      vgeo(~idnan,j,t) = g./f(~idnan,j).*nansum(Dxnan.*sshj,2);

    end

  end

end

end % end main function


function [D] = fdmatrix(x,diff_ord,accu_ord)

% Function to output a finite difference matrix D based on the function fdcoefs
%   x is the grid of nodes [x0 x1 x2...xN]
%   diff_ord is the differentiation order
%   accu_ord is the accuracy order

n = accu_ord+diff_ord-1;    %n is the number of points in the stencil
m = diff_ord;               %m is the differentiation order

npoints = length(x);        %Number of grid points
nedge = floor((n+1)/2);     %Number of end points requiring consideration 

D = zeros(npoints);

if mod(n,2)==1
    for i = 1:nedge
        D(i,1:(n+1)) = fdcoefs(m,n,x(1:(n+1)),x(i));
    end
    for i = (nedge+1):(npoints-nedge)
        D(i,(i-nedge):(i+nedge-1)) = fdcoefs(m,n,x((i-nedge):(i+nedge-1)),x(i));
    end
    for i = (npoints-nedge+1):(npoints)
        D(i,(end-n):end) = fdcoefs(m,n,x((end-n):end),x(i));
    end
else
    for i = 1:nedge
        D(i,1:(n+1)) = fdcoefs(m,n,x(1:(n+1)),x(i));
    end
    for i = (nedge+1):(npoints-nedge)
        D(i,(i-nedge):(i+nedge)) = fdcoefs(m,n,x((i-nedge):(i+nedge)),x(i));
    end
    for i = (npoints-nedge+1):(npoints)
        D(i,(end-n):end) = fdcoefs(m,n,x((end-n):end),x(i));
    end
end

%D = sparse(D);

end % end functioon fdmatrix


function coefs = fdcoefs(m,n,x,xi)

% Finite difference weights (Fornberg)
%   m: Differentiation order
%   n: size of the stencil (3,5,7,9)
%   x: coordinates of the points in the stencil
%   xi: coordinate of the evaluation point
% Number of points in the formula n+1: formal order n-m+1 (irregular points)

c1= 1;
c4= x(1)-xi;

c= zeros(n+1,m+1);
c(1,1)= 1;

for i=1:n-1;
    mn= min([i,m]);
    c2= 1;
    c5= c4;
    c4= x(i+1)-xi;
    for j=0:i-1;
        c3= x(i+1)-x(j+1);
        c2= c2*c3;
        for k= mn:-1:1;
            c(i+1,k+1)= c1*(k*c(i,k)-c5*c(i,k+1))/c2;
        end;
        c(i+1,1)= -c1*c5*c(i,1)/c2;
        for k=mn:-1:1;
            c(j+1,k+1)= (c4*c(j+1,k+1)-k*c(j+1,k))/c3;
        end;
        c(j+1,1)= c4*c(j+1,1)/c3;
    end;
    c1= c2;
end;

coefs= c(:,m+1)';

end % end function fdcoefs

