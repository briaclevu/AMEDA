% plot_tracking.m
%
% Make a movie with tracks
%
%=========================

start
clear; clc;

movie=0;
h_coast=1;
bath=1;
track=0; % add satellite track (1)
temp=0;% use sst AVHRR
chlo=0;% use chlo composite
buoy=0;argo=1;glider=0; % add buoy,argo,gldier trajectories

%----------------------------------------
% source of data driving the netcdf format
source = 'AVISO';

%----------------------------------------
% configuration to be used has reference of the v1 of the atlas
config = 'DYNED_MED_cyclo';
year='2000_2017';
name={'LEV','BWE'};

% Definition of the parameters specific for the experiment is needed 
%mod_eddy_params(['keys_sources_',source,'_',config])
run(['keys_sources_',source,'_',config])
load('param_eddy_tracking')
path_in2='/home/blevu/DATA/AVISO/DYNED_MED/';
nc_ssh=[path_in2,'ssh_adt_DYNED_MED_2000_2017.nc'];

%% plot a movie or pdf
% make figure's dirs
system(['mkdir ',path_out,'figures']);

% load bathymetry
if bath
    load('bathy/bathy_med')
end

% load high coastal def
if h_coast
    coast=csvread('bathy/new_bathy.csv');
end

% load buoy trajectory
if buoy
    dirdata=[path_data,'drifters/'];
    list=dir([dirdata,'*.mat']);

    i=1;
    for j=1:length(list)
        file = [dirdata,list(j).name];
        if exist(file,'file')
            load(file)
            if ~isempty(dateb)
                drifter(i).name=list(j).name(1:end-4);
                drifter(i).date=dateb;
                drifter(i).lat=latb;
                drifter(i).lon=lonb;
                i=i+1;
            end
        end
    end
    % duration for the ruban
    nbb=15;
end

if argo
    input_dir = '/home/blevu/DATA/CORIOLIS/argo_profiles/med/argo_profiles_Med_all_modes.nc';
    input_dir = '/home/blevu/DATA/CORIOLIS/argo_profiles/med/coloc_anoms_Med_all_modes.nc';
    time_1_argo = '01-Jan-2000';
    time_1_track = '01-Jan-2000';
    time_argo = floor(double(ncread(input_dir,'time_step')+datenum(time_1_argo)-1));
    %ID_traj = double(ncread(input_dir,'ID_track'));
    ID_float = double(ncread(input_dir,'ID_float'));
    %IsInEddy = double(ncread(input_dir,'IsInEddy'));
    Xargo = double(ncread(input_dir,'Xargo'));
    Yargo = double(ncread(input_dir,'Yargo'));

    %
    dirdata=[path_data,'argo_trajs/'];
    list=dir([dirdata,'*.mat']);
    i=1;
    for j=1:length(list)
        file = [dirdata,list(j).name];
        if exist(file,'file')
            load(file)
            if ~isempty(dateb)
                argos(i).name=list(j).name(17:end-4);
                argos(i).date=dateb;
                argos(i).lat=latb;
                argos(i).lon=lonb;
                i=i+1;
            end
        end
    end
    % duration for the ruban
    nba=20;
end

% load tracking result
load([path_out,'eddy_tracks2_',year])  
tracks=tracks2;

% load grid
eval(['[x,y,mask,~,~,~] = load_fields_',source,'(1,1,deg);'])

% load SMDT
if strcmp(config,'PRO2016_sla') || strcmp(config,'PRO2017_sla')
    mdt=ncread([path_in,'SMDT-MED-2014-REF20.nc'],'mdt')';
    lon_mdt=ncread([path_in,'SMDT-MED-2014-REF20.nc'],'lon')';
    lat_mdt=ncread([path_in,'SMDT-MED-2014-REF20.nc'],'lat')';
    mdt1=interp2(lon_mdt,lat_mdt,mdt,x,y);
end

% set area
for named=name
minlat=min(y(:));
maxlat=max(y(:));
minlon=min(x(:));
maxlon=max(x(:));
cmin=-25;
cmax=25;
sp=3;
coeff=2;
% set boundaries
if strcmp(named{1},'BWE')
    minlat=35;
    maxlat=45;
    minlon=-6;
    maxlon=16;
    minx=330;
    maxx=700;
    miny=330;
    maxy=730;
    tmax=[26 21 18 18 20];
    tmin=[20 15 15 12 12 12];
    cmin=-25;
    cmax=25;
elseif strcmp(named{1},'LEV')
    minlat=30;
    maxlat=40;
    minlon=14;
    maxlon=36;
    cmin=-25;
    cmax=25;
elseif strcmp(named{1},'ALG')
    minlat=36;
    maxlat=41;
    minlon=-.5;
    maxlon=9.5;
    cmin=-25;
    cmax=25;
elseif strcmp(named{1},'LYB')
    minlat=31;
    maxlat=36;
    minlon=20;
    maxlon=30;
    cmin=-25;
    cmax=25;
end
system(['mkdir ',path_out,'/figures/',named{1}])
system(['mkdir ',path_out,'/figures/',named{1},'_ID'])

% set periode
dayi=[2000 01 01 12 0 0];
%dstp=datenum(dayi)-datenum([2000 01 01 12 0 0]);

if temp
    stp_list = [610 1 1]; % temp
    stp_list = [400 1 1]; % temp
elseif movie
    stp_list =[1 367 732 1097 1462 1828 2193 2558 2923 3289 3654 4019 4384 4750 5115 5480 stepF+1];
else
    stp_list =[stepF-12 stepF];    
    stp_list =[1 366 731 1097 stepF];
    stp_list =[1 367 732 1097 1462 1828 2193 2558 2923 3289 3654 4019 4384 4750 5115 5480 5845 6211 stepF+1];
    fig_list =[640 666 815 849 854 881];
end
dstp=3653;

% set the limit of eddies life
cut=0;

m_proj('Mercator','lat',[minlat maxlat],'lon',[minlon maxlon]);

close all

for n=1:length(stp_list)-1
    system(['mkdir ',path_out,'/figures/',named{1},'/',num2str(n+1999)])
    system(['mkdir ',path_out,'/figures/',named{1},'_ID/',num2str(n+1999)])
    %dstp = stp_list(n)-1; 
    dstp = 0;
    
    if movie
    % prepare fig
        hfig=figure('visible','off');
        set(hfig,'Position',[0 0 1000 600])
        if temp
            M=VideoWriter([path_out,'figures/sst_',named{1},'_t',num2str(stp_list(n)),'_t',num2str(stp_list(n+1)-1),'_article.avi']);
        elseif strcmp(config,'PRO2017_adt')
            M=VideoWriter([path_out,'figures/adt_',config,'_',named{1},'_',year,'.avi']);
        elseif strcmp(config,'PRO2017_sla')
            M=VideoWriter([path_out,'figures/sla_',config,'_',named{1},'_',year,'.avi']);
        end
        M.FrameRate = 1;
        open(M);
    end
    % load avhrr sst
    if temp
        if n==1
            sstfile=['/home/blevu/DATA/PODAAC/SST_MED_2015_121_365.nc'];
        elseif n==2
            sstfile=['/home/blevu/DATA/PODAAC/SST_MED_2016_001_366.nc'];
        elseif n==3
            sstfile=['/home/blevu/DATA/PODAAC/SST_MED_2017_001_120.nc'];
        end
        sst=ncread(sstfile,'sea_surface_temperature')-273.15;
        sst_lon=ncread(sstfile,'lon');
        sst_lon=sst_lon(minx:maxx,miny:maxy);
        sst_lat=ncread(sstfile,'lat');
        sst_lat=sst_lat(minx:maxx,miny:maxy);
        sst_time=datevec(ncread(sstfile,'time')/60/60/24+...
            datenum('1981-01-01 00:00:00'));
        %
        sst(sst<=12)=nan;
    end
    
    for day=stp_list(n):stp_list(n+1)-1%stp_list(n):stp_list(n+1)-1%1707:1712%[1291:1389 2315:3363 3713:3799 4038:4081 4174:4297 4787:5579]%[1180:1796,2565:2670,3420:5091]%%98:104%[1 24 43 151 212 251 271]
        
        close all
        if ~movie
            figure, close
            hfig=figure('visible','off');
            set(hfig,'Position',[0 0 1000 600])
            if strcmp(name,'Indi')
                set(hfig,'Position',[0 0 1000 600])
            elseif strcmp(name,'Oman')
                set(hfig,'Position',[0 0 1000 400])
            elseif strcmp(name,'Aden')
                set(hfig,'Position',[0 0 900 600])
            elseif strcmp(name,'West')
                set(hfig,'Position',[0 0 1000 900])
            end
        end
    % reference day at 12:00
        dr = datevec(datenum(dayi)+day-1);
    % pcolor with wwh
        eval(['[x,y,mask,u1,v1,~] = load_fields_',source,'(day+dstp,1,deg);'])
        eval(['[x1,y1,mask1,~,~,~] = load_fields_',source,'(day+dstp,1,3);'])
        ssh1 = squeeze(permute(ncread(nc_ssh,s_name,[1 1 day+dstp],[Inf Inf 1]),[2,1,3]));
        xnan=x;
        ynan=y;
        vel=sqrt(u1.^2+v1.^2);
        xnan(isnan(vel) | vel<0.05) = nan;
        ynan(isnan(vel) | vel<0.05) = nan;
        m_contourf(x1,y1,ssh1*100,50)
        shading flat
        h = colorbar;
        caxis([cmin cmax])
        colormap(soft(40,0.2))
        %colormap(cbrewer('seq','Oranges',20))
        % add quiver velocities
        hold on
        m_quiver(xnan(1:sp:end,1:sp:end),ynan(1:sp:end,1:sp:end),...
            u1(1:sp:end,1:sp:end),v1(1:sp:end,1:sp:end),coeff,'color',[0.3 0.3 0.3])
        % put a land mask
        if h_coast
            [X,~]=m_ll2xy(coast(:,1),coast(:,2),'clip','patch');
            k=find(isnan(X(:,1)));
            for i=1:length(k)-1,
                xc=coast([k(i)+1:(k(i+1)-1) k(i)+1],1);
                yc=coast([k(i)+1:(k(i+1)-1) k(i)+1],2);
                m_patch(xc,yc,[.9 .9 .9]);
            end
        else
            m_coast('patch',[.9 .9 .9],'edgecolor','k');
        end
        % grid and fancy the map
        m_grid('tickdir','in','linewidth',1,'linestyle','none');
        set(gcf,'color','w'); % set figure background to white
        
        % add buoy trajectory
        if buoy
            for i=1:length(drifter)
                dateb=datenum(drifter(i).date)-datenum(dr);
                latb=drifter(i).lat;
                lonb=drifter(i).lon;
                ruban=find(dateb<0 & dateb>-nbb &...
                    lonb<maxlon & lonb>minlon & latb<maxlat & latb>minlat);
                if ~isempty(ruban)
                    L=length(ruban);
                    m_scatter(lonb(ruban),latb(ruban),10,...
                        [ones(L,1) -dateb(ruban)/nbb -dateb(ruban)/nbb],...
                        'filled')
                end
                ruban=find(dateb>0 & dateb<nbb &...
                    lonb<maxlon & lonb>minlon & latb<maxlat & latb>minlat);
                ruban=[];% no futur trajectories
                if ~isempty(ruban)
                    L=length(ruban);
                    m_scatter(lonb(ruban),latb(ruban),10,...
                        [1-dateb(ruban)/nbb zeros(L,1) zeros(L,1)],...
                        'filled')
                end
                mrub=find(abs(dateb)==min(abs(dateb)),1);
                if abs(dateb(mrub))<=1 & abs(dateb(mrub))>=-1 
                    m_plot(lonb(mrub),latb(mrub),'ok','markerfacecolor',[1 0 0],...
                        'markersize',5)
                end
            end
            %m_text(maxlon-2.4,maxlat-0.65,['DRIFTER + and - ',num2str(nbb),' days (red)'],'FontSize',8)
            m_text(maxlon-5,maxlat+0.6,['DRIFTERS last ',num2str(nbb),' days (red circles)'],'FontSize',6)
        end
        
        
        % satellite tracks of the day (use get_tracks.sh in DATA/nrt)
        % (time,lat,lon,flag,track,cycle)
        if track
            nsat=0;
            nday=2;
            for sat={'al','h2','c2','j2','j3','e1','e2','en','enn','g2','j1','j1g','j1n','tp','tpn'}
                for i=day-nday:day+nday
                    dri=datevec(datenum(dayi)+i-1);
                    trackname=[path_tracks1,char(sat),'/dt_med_',char(sat),...
                            '_adt_vfec_',datestr(dri,'yyyymmdd'),'.nc'];
                    if ~exist(trackname,'file')
                        trackname=[path_tracks2,char(sat),'/nrt_med_',char(sat),...
                                '_adt_vfec_',datestr(dri,'yyyymmdd'),'.nc'];
                    end
                    if exist(trackname,'file')
                        lont=ncread(trackname,'longitude');
                        lont(lont>180)=lont(lont>180)-360;
                        latt=ncread(trackname,'latitude');
                        if ~isempty(lont(lont>minlon & lont<maxlon & latt>minlat & latt<maxlat))
                            nsat=nsat+1;
                            if strcmp(sat,'h2')
                                m_plot(lont,latt,'xk','markersize',4)                                
                            else
                                m_plot(lont,latt,'xk','markersize',4)
                            end
                        end
                    end
                end
            end
            m_text(maxlon-5,maxlat+0.2,[num2str(nsat),' tracks +/- ',num2str(nday),' days (cross)'],'FontSize',6)
        end
        
        % add eddy tracking (use get_aviso_nrt.sh in DATA/nrt -> mk_nc_along.m)
        for i=1:length(tracks)
            ind=[];
            ind=find(tracks(i).step==day,1);
            if ~isempty(ind)
                %CD = [(tracks(i).x1)';(tracks(i).y1)'];
                CD = [(tracks(i).xbary1)';(tracks(i).ybary1)'];
                dura = tracks(i).step(ind)-tracks(i).step(1)+1;
                if tracks(i).step(end)-tracks(i).step(1)+1>=cut
                    if dura<=30
                        m_plot(CD(1,1:ind),CD(2,1:ind),'-','color',[0.4 0.4 0.4],'linewidth',1.5)
                        m_plot(CD(1,1),CD(2,1),'sk','MarkerFaceColor','k','MarkerSize',4)
                    else
                        m_plot(CD(1,max(1,ind-30):ind),CD(2,max(1,ind-30):ind),'-','color',[0.4 0.4 0.4],'linewidth',1.5)
                        m_plot(CD(1,max(1,ind-30)),CD(2,max(1,ind-30)),'sk','MarkerFaceColor',[0.3 0.3 0.3],'MarkerSize',4)
                    end
                    m_plot(CD(1,ind),CD(2,ind),'ok','MarkerFaceColor','k',...
                        'MarkerSize',min(10,round(dura/90+4)))
                    %
                    lonlat3=tracks(i).shapes3{ind};
                    m_plot(lonlat3(1,:),lonlat3(2,:),'--k','linewidth',1)
                    %
                    lonlat1=tracks(i).shapes1{ind};
                    if ~isnan(lonlat1)
                        if tracks(i).type(1)==-1
                            col=[.1 .1 .9]; % anticyclone
                        else
                            col=[.9 .1 .1]; % cyclone
                        end
                        if tracks(i).large1(ind)==0
                            wid='-'; % largest eddy
                        else
                            wid='--'; % true eddy
                        end
                        m_plot(lonlat1(1,:),lonlat1(2,:),wid,'color',col,'linewidth',1.5)
                    end
                    %
                    lonlat2=tracks(i).shapes2{ind};
                    if ~isnan(lonlat2)
                        list=[tracks(i).interaction(ind) tracks(i).interaction2(ind)];
                        list=list(~isnan(list));
                        for j=list
                            indj = find(tracks(j).step==day+dstp,1);
                            if ~isempty(indj)
                                col= [.1 .5 .1]; % interaction
                                if tracks(i).split(ind)==1 || tracks(i).merge(ind)==1 ||...
                                    tracks(i).split2(ind)==1 || tracks(i).merge2(ind)==1 ||...
                                    tracks(j).split(indj)==1 || tracks(j).merge(indj)==1 ||...
                                    tracks(j).split2(indj)==1 || tracks(j).merge2(indj)==1
                                    if (tracks(i).split(ind)==1 && tracks(i).merge(ind)==1) ||...
                                        (tracks(i).split2(ind)==1 && tracks(i).merge2(ind)==1) ||...
                                        (tracks(j).split(indj)==1 && tracks(j).merge(indj)==1) ||...
                                        (tracks(j).split2(indj)==1 && tracks(j).merge2(indj)==1)
                                            col= [1 1 1]; % merging + splitting
                                    else
                                            col= [.1 .9 .1]; % merging or splitting
                                    end
                                end
                                m_plot(lonlat2(1,:),lonlat2(2,:),'-','color',col,'linewidth',1) % double eddy
                            end
                        end
                    end
                    %
                    if CD(1,ind)<maxlon && CD(1,ind)>minlon && CD(2,ind)>minlat && CD(2,ind)<maxlat
                        m_text(CD(1,ind),CD(2,ind)-.1,['  ',num2str(dura)],...
                            'color',[0 0 0],'FontSize',min(10,round(dura/90+4)),'fontWeight','bold')
                        %m_text(CD(1,ind),CD(2,ind)-.1,[num2str(i)],...
                         %   'color',[.5 .5 .5],'FontSize',8,'fontWeight','bold')
                    end
                else
                    m_plot(CD(1,ind),CD(2,ind),'ok',...
                        'MarkerFaceColor','w','MarkerSize',4)
                end
            end
        end
        
        % add argo trajectory
        if argo
%             IND = find(time_argo==floor(datenum(dr)));
%             if ~isempty(IND)
%                 for i=1:length(IND)
%                     if IsInEddy(IND(i))==1 || IsInEddy(IND(i))==-1
%                         m_plot(Xargo(IND(i)),Yargo(IND(i)),'dk','markerfacecolor',[1 1 0],...
%                         'markersize',8)
%                     end
%                 end
%             end
            for i=1:length(argos)
                dateb=datenum(argos(i).date)-datenum(dr);
                lonb=argos(i).lon;
                latb=argos(i).lat;
                mrub=find(dateb(dateb<0)==max(dateb(dateb<0)));
                ruban=find(dateb<0 & dateb>-nba &...
                    diff([dateb;10000]) > 0.5 &...
                    lonb<maxlon & lonb>minlon & latb<maxlat & latb>minlat);
                if ~isempty(ruban)
                    L=[length(ruban(:)) i];
                    m_plot(lonb(ruban(:)),latb(ruban(:)),'-','color',[1 0 1],'linewidth',1)
                    m_plot(lonb(ruban(:)),latb(ruban(:)),'dk','markerfacecolor',[1 0 1],...
                        'markersize',4)
                    %m_text(lonb(ruban(1)),latb(ruban(1))+0.1,argos(i).name)
                end
                if dateb(mrub)>-1
                    m_plot(lonb(mrub(end)),latb(mrub(end)),'dk','markerfacecolor',[1 0 1],...
                        'markersize',6)

                end
            end
            %m_text(maxlon-2.8,maxlat-0.35,['ARGO + and - ',num2str(nba),' days (purple)'],'FontSize',8)
            %m_text(maxlon-5,maxlat+0.4,['ARGO last ',num2str(nba),' days (pink diamonds)'],'FontSize',6)
        end
        % add legend on the map and information
        %set(get(h,'title'),'string','ADT (cm)','FontSize',8,'fontWeight','bold')
        
        hold off
        % add legend on the map and information
        if temp
            set(get(h,'title'),'string','T(degC)','FontSize',8,'fontWeight','bold')
            if sst_time(stp,4)==2 || sst_time(stp,4)==13
                m_text(maxlon-2,minlat+.2,[datestr(sst_time(stp,:),'dd-mmm-yyyy HH:MM'),' | NPP'],'FontSize',10,'fontWeight','bold');
            else
                m_text(maxlon-2,minlat+.2,[datestr(sst_time(stp,:),'dd-mmm-yyyy HH:MM'),' | METOP'],'FontSize',10,'fontWeight','bold');
            end
        else
            if strcmp(config,'PRO2017_adt')
                set(get(h,'title'),'string','ADT (cm)','FontSize',8,'fontWeight','bold')
            elseif strcmp(config,'PRO2017_sla')
               set(get(h,'title'),'string','SLA (cm)','FontSize',8,'fontWeight','bold')
            end
            set(get(h,'title'),'string','ADT (cm)','FontSize',8,'fontWeight','bold')
            m_text(minlon+.2,maxlat+.2,[datestr(dr,1)],...%' / degradation = ',num2str(deg)],...
                'FontSize',10,'fontWeight','bold')
        end

        %prepare the print in a pdf
        if movie
            figinfo = hardcopy(hfig,'-dzbuffer','-r0');
            writeVideo(M,im2frame(figinfo));
            %frame = getframe(gcf);
            %writeVideo(M,frame);
        else
            set(hfig,'PaperPosition',[-1.3,-0.4,12,7])
            set(hfig,'PaperSize',[9.5,8.5])
            set(hfig,'PaperPosition',[-1.3,-0.4,10,7])
            set(hfig,'PaperSize',[7.5,6.5])
            %print(hfig,[path_out,'figures/',named{1},'/',num2str(n+1999),'/tracking_',config,'_',named{1},'_',datestr(dr,'yyyymmdd')],'-dpng','-r200')
            print(hfig,[path_out,'figures/',named{1},'/',num2str(n+1999),'/tracking_',config,'_',named{1},'_',datestr(dr,'yyyymmdd')],'-dpng','-r300')
        end
        
        % figure with ID
        if ~movie
            figure, close
            hfig=figure('visible','off');
            set(hfig,'Position',[0 0 1000 600])
            if strcmp(name,'Indi')
                set(hfig,'Position',[0 0 1000 600])
            elseif strcmp(name,'Oman')
                set(hfig,'Position',[0 0 1000 400])
            elseif strcmp(name,'Aden')
                set(hfig,'Position',[0 0 900 600])
            elseif strcmp(name,'West')
                set(hfig,'Position',[0 0 1000 900])
            end
        end
        % reference day at 12:00
        dr = datevec(datenum(dayi)+day-1);
        % put a land mask
        m_coast('color','k');
        % grid and fancy the map
        m_grid('tickdir','in','linewidth',1,'linestyle','none');
        %colormap(cbrewer('seq','Oranges',20))
        set(gcf,'color','w'); % set figure background to white
        hold on

        % add eddy center (use get_aviso_nrt.sh in DATA/nrt -> mk_nc_along.m)
        for i=1:length(tracks)
            ind=[];
            ind=find(tracks(i).step==day+dstp,1);
            if ~isempty(ind)
                CD = [(tracks(i).xbary1)';(tracks(i).ybary1)'];
                dura = tracks(i).step(ind)-tracks(i).step(1)+1;
                if tracks(i).step(end)-tracks(i).step(1)+1>=cut
                    %
                    lonlat1=tracks(i).shapes1{ind};
                    if ~isnan(lonlat1)
                        if tracks(i).type(1)==-1
                            col=[.1 .1 .9]; % anticyclone
                        else
                            col=[.9 .1 .1]; % cyclone
                        end
                    end
                    %
                    if CD(1,ind)<maxlon && CD(1,ind)>minlon && CD(2,ind)>minlat && CD(2,ind)<maxlat
                        %m_text(CD(1,ind),CD(2,ind)-.1,['  ',num2str(dura)],...
                         %   'color',[0 0 0],'FontSize',min(10,round(dura/90+4)),'fontWeight','bold')
                        m_text(CD(1,ind),CD(2,ind)-.1,['  ',num2str(i)],...
                            'color',col,'FontSize',min(10,round(dura/90+4)),'fontWeight','bold')
                        m_plot(CD(1,ind),CD(2,ind),'ok','MarkerFaceColor',col,...
                            'MarkerSize',min(10,round(dura/90+4)))
                    end
                end
            end
        end
        
        % add argo name
        if argo
            for i=1:length(argos)
                dateb=datenum(argos(i).date)-datenum(dr);
                lonb=argos(i).lon;
                latb=argos(i).lat;
                mrub=find(dateb(dateb<0)==max(dateb(dateb<0)));
                ruban=find(dateb<0 & dateb>-nba &...
                    diff([dateb;10000]) > 0.5 &...
                    lonb<maxlon & lonb>minlon & latb<maxlat & latb>minlat);
                if ~isempty(ruban)
                    m_plot(lonb(ruban(1)),latb(ruban(1)),'dk','markerfacecolor',[1 0 1],...
                        'markersize',6)
                    m_text(lonb(ruban(1)),latb(ruban(1))+0.15,argos(i).name,...
                        'color',[1 0 1],'FontSize',6,'fontWeight','bold')
                end
            end
        end

        % add legend on the map and information
        hold off
        % add legend on the map and information
        m_text(minlon+.2,maxlat+.2,[datestr(dr,1)],...%' / degradation = ',num2str(deg)],...
            'FontSize',10,'fontWeight','bold')

        %prepare the print in a pdf
        set(hfig,'PaperPosition',[-1.3,-0.4,10,7])
        set(hfig,'PaperSize',[7.5,6.5])
        print(hfig,[path_out,'figures/',named{1},'_ID/',num2str(n+1999),'/tracking_ID_',config,'_',named{1},'_',datestr(dr,'yyyymmdd')],'-dpng','-r150')

    end
    
    if movie, close(M); end
        
end

end
