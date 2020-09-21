% Compute geostrophic velocities from ssh serie
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User modifications
sshname='croco'; % choose adt or sla
domname='CRE';
postname='_20150101_20150110';

% directory
dirin=[pwd,'/',domname,'/'];
mkdir(dirin,'geo')

n=9; % n-stencil finite difference centred method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% read fields
ssh=ncread([dirin,'ssh_',sshname,'_',domname,postname,'.nc'],'ssh');
[N,M,L]=size(ssh);

lon=ncread([dirin,'lon_lat_',sshname,'_',domname,'.nc'],'lon');
lat=ncread([dirin,'lon_lat_',sshname,'_',domname,'.nc'],'lat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute geostrophic velocities stencil-9

[ugeo,vgeo] = compute_speed_geo(n,lon,lat,ssh);

nccreate([dirin,'geo/ssu_',sshname,'_',domname,'_geo',postname,'.nc'],'day','Dimension',{'day' Inf});
nccreate([dirin,'geo/ssu_',sshname,'_',domname,'_geo',postname,'.nc'],'ugeo','Dimension',{'x' N 'y' M 'day' L});
ncwrite([dirin,'geo/ssu_',sshname,'_',domname,'_geo',postname,'.nc'],'ugeo',ugeo);

nccreate([dirin,'geo/ssv_',sshname,'_',domname,'_geo',postname,'.nc'],'day','Dimension',{'day' Inf});
nccreate([dirin,'geo/ssv_',sshname,'_',domname,'_geo',postname,'.nc'],'vgeo','Dimension',{'x' N 'y' M 'day' L});
ncwrite([dirin,'geo/ssv_',sshname,'_',domname,'_geo',postname,'.nc'],'vgeo',vgeo);




