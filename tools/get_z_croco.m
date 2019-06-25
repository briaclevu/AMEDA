function z = get_z_croco(ncname,nc_ssh,tindex,theta_s,theta_b,hc,N,type,dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function z = get_z_croco(ncname,nc_ssh,tindex,hc,N,type,dim);
%
%  this function compute the depth of rho or w points for CROCO
%  in case of the new v transform (Shcheptekin, 2006)
%
%  On Input:
%
%	 tindex: step index
%    hc:     critical depth
%    N:      number of level (0 to read it from saved S-curve)
%    type:   'r': rho point 'w': w point
%    dim:	 dimension arrangement
%
%  On Output:
%
%    z       Depths (m) of RHO- or W-points (3D matrix).
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read bathymetry from grid
h=squeeze(ncread(ncname,'h'));
[M,L]=size(h);

% read zeta in hourly 3D output
zeta=mean(ncread(nc_ssh,'ssh',[1,1,max(1,(tindex-1)*24)],[Inf,Inf,24]),3);

%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.

if  ~N
% read S-Curves and the number of level

   if type=='w'
      sc_w=squeeze(ncread(ncname,'sc_w'));
      Cs_w=squeeze(ncread(ncname,'Cs_w'));
      N=length(sc_w);
   else
      sc_r=squeeze(ncread(ncname,'sc_r'));
      Cs_r=squeeze(ncread(ncname,'Cs_r'));
      N=length(sc_r);
  end

else
% Compute S-curves

   sc_r=zeros(N,1);
   Cs_r=zeros(N,1);
   sc_w=zeros(N+1,1);
   Cs_w=zeros(N+1,1);
   %
   ds=1./N;

   if type=='w'
      sc_w(1) = -1.0;
      sc_w(N+1) =  0;
      Cs_w(1) = -1.0;
      Cs_w(N+1) =  0;
      %
      sc_w(2:N) = ds*([1:N-1]-N);
      Cs_w=csf(sc_w, theta_s,theta_b);
      N=N+1;
   else
      sc= ds*([1:N]-N-0.5);    
      Cs_r=csf(sc, theta_s,theta_b);
      sc_r=sc;
   end

end

%
% Create S-coordinate system: based on model topography h(i,j),
% fast-time-averaged free-surface field and vertical coordinate
% transformation metrics compute evolving depths of of the three-
% dimensional model grid. Also adjust zeta for dry cells.

h(h==0)=1.e-2;
Dcrit=0.01;   % min water depth in dry cells
zeta(zeta<(Dcrit-h))=Dcrit-h(zeta<(Dcrit-h));

hinv=1./h;
z1=zeros(N,M,L);

if type=='w'
    cff1=Cs_w;
    cff2=sc_w+1;
    sc=sc_w;
else
    cff1=Cs_r;
    cff2=sc_r+1;
    sc=sc_r;
end

h2=(h+hc);
cff=hc*sc;
h2inv=1./h2;
for k=1:N
    z0=cff(k)+cff1(k)*h;
    z1(k,:,:)=z0.*h./(h2) + zeta.*(1.+z0.*h2inv);
end

z = permute(z1,dim);

return

            
        

