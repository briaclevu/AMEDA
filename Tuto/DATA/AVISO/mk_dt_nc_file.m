% Arrange and rename mapped data field coming from CMEMS website 
% to be computed in mod_eddy_* routines

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User modifications
domname='MED';

dirin=[pwd,'/phy_l4/'];

dirout=[pwd,'/',domname,'/'];

system(['mkdir ',dirout])

% data from aviso.altimetry.fr
inname='dt_med_allsat_phy_l4';

M=textread([pwd,'/date.txt']);
datename=['_',num2str(M(1)),'_',num2str(M(end))];

% build outfile name
outname=['dt_',domname,datename];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the fields (lon,lat,time, ssh, u and v)
ssh=ncread([dirin,inname,datename,'.nc'],'adt');
u=ncread([dirin,inname,datename,'.nc'],'ugos');
v=ncread([dirin,inname,datename,'.nc'],'vgos');

time=ncread([dirin,inname,datename,'.nc'],'time');
lon=ncread([dirin,inname,datename,'.nc'],'longitude');
lat=ncread([dirin,inname,datename,'.nc'],'latitude');

% area med
lon1=-6-360;
lon2=37-360;
lat1=30;
lat2=46;

% Resize for the selected area
u1=u((lon>=lon1+360 & lon<=lon2+360),(lat>=lat1 & lat<=lat2)',:);

[N1,M1,L1]=size(u1);
mask1=squeeze(u1(:,:,1)*0+1);
mask1(isnan(mask1))=0;

% Arrange longitude (360 -> OE)
lon2d=double(repmat(lon(lon>=lon1+360 & lon<=lon2+360),[1 M1]));
lat2d=double(repmat(lat(lat>=lat1 & lat<=lat2)',[N1 1]));

% Create grid file
if ~exist([dirout,'lon_lat_',domname,'.nc'],'file')
  nccreate([dirout,'lon_lat_',domname,'.nc'],'lon','Dimensions',{'x' N1 'y' M1});
  nccreate([dirout,'lon_lat_',domname,'.nc'],'lat','Dimensions',{'x' N1 'y' M1});
  nccreate([dirout,'lon_lat_',domname,'.nc'],'mask','Dimensions',{'x' N1 'y' M1});
  ncwrite([dirout,'lon_lat_',domname,'.nc'],'mask',mask1);
  ncwrite([dirout,'lon_lat_',domname,'.nc'],'lon',lon2d);
  ncwrite([dirout,'lon_lat_',domname,'.nc'],'lat',lat2d);
end

% Create fields file in the right box
system(['ncks -O ',dirin,inname,datename,'.nc -v longitude,latitude,time,adt ',...
    '-d longitude,',num2str(lon1+360),'.0,',num2str(lon2+360),...
    '.0 -d latitude,',num2str(lat1),'.0,',num2str(lat2),'.0 ',...
    dirout,'ssh_',outname,'.nc']);
system(['ncks -O ',dirin,inname,datename,'.nc -v longitude,latitude,time,ugos ',...
    '-d longitude,',num2str(lon1+360),'.0,',num2str(lon2+360),...
    '.0 -d latitude,',num2str(lat1),'.0,',num2str(lat2),'.0 ',...
    dirout,'ssu_',outname,'.nc']);
system(['ncks -O ',dirin,inname,datename,'.nc -v longitude,latitude,time,vgos ',...
    '-d longitude,',num2str(lon1+360),'.0,',num2str(lon2+360),...
    '.0 -d latitude,',num2str(lat1),'.0,',num2str(lat2),'.0 ',...
    dirout,'ssv_',outname,'.nc']);

% Change name to be used by the mod_eddy routines
system(['ncrename -v time,day -v adt,ssh ',dirout,'ssh_',outname,'.nc']);
system(['ncrename -v time,day -v ugos,u ',dirout,'ssu_',outname,'.nc']);
system(['ncrename -v time,day -v vgos,v ',dirout,'ssv_',outname,'.nc']);
