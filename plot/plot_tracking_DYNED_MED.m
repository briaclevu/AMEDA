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
export_out=0;% export results

%----------------------------------------
% source of data driving the netcdf format
source = 'AVISO';

%----------------------------------------
% configuration to be used has reference of the v1 of the atlas
config = 'DYNED_MED_cyclo';
year='2000_2017';
name={'BWE'};

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
    list=unique(ID_float);
    i=1;
    clear argos
    for j=1:length(list)
        IND=find(ID_float==list(j));
        argos(i).name = list(j);
        argos(i).date = time_argo(IND);
        argos(i).lon =Xargo(IND);
        argos(i).lat = Yargo(IND);
        [argos(i).date,I] = sort(argos(i).date);
        argos(i).lon = argos(i).lon(I);
        argos(i).lat = argos(i).lat(I);
        i=i+1;
    end

    %
%     dirdata=[path_data,'argo_trajs/'];
%     list=dir([dirdata,'*.mat']);
%     i=1;
%     for j=1:length(list)
%         file = [dirdata,list(j).name];
%         if exist(file,'file')
%             load(file)
%             if ~isempty(dateb)
%                 argos(i).name=list(j).name(17:end-4);
%                 argos(i).date=dateb;
%                 argos(i).lat=latb;
%                 argos(i).lon=lonb;
%                 i=i+1;
%             end
%         end
%     end
    % duration for the ruban
    nba=15;
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
mkdir([path_out,'/figures'],named{1})
mkdir([path_out,'/figures'],[named{1},'_ID'])

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

for n=11:length(stp_list)-1
    mkdir([path_out,'/figures/',named{1}],num2str(n+1999))
    mkdir([path_out,'/figures/',named{1},'_ID'],num2str(n+1999))
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
        % set the scale
        if strcmp(named{1},'BWE')
            u1(331,280)=1;
            v1(331,280)=0;
        elseif strcmp(named{1},'LEV')
            u1(169,850)=1;
            v1(169,850)=0;
        end
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
            u1(1:sp:end,1:sp:end)/coeff,v1(1:sp:end,1:sp:end)/coeff,...
            'autoscale','off','color',[.3 .3 .3])
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
                        m_plot(CD(1,1),CD(2,1),'sk','MarkerFaceColor','k','MarkerSize',3)
                    else
                        m_plot(CD(1,max(1,ind-30):ind),CD(2,max(1,ind-30):ind),'-','color',[0.4 0.4 0.4],'linewidth',1.5)
                        m_plot(CD(1,max(1,ind-30)),CD(2,max(1,ind-30)),'sk','MarkerFaceColor',[0.3 0.3 0.3],'MarkerSize',3)
                    end
                    m_plot(CD(1,ind),CD(2,ind),'ok','MarkerFaceColor','k',...
                        'MarkerSize',min(8,round(dura/120+3)))
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
                            'color',[0 0 0],'FontSize',min(8,round(dura/120+4)),'fontWeight','bold')
                        %m_text(CD(1,ind),CD(2,ind)-.1,[num2str(i)],...
                         %   'color',[.5 .5 .5],'FontSize',8,'fontWeight','bold')
                    end
                else
                    m_plot(CD(1,ind),CD(2,ind),'ok',...
                        'MarkerFaceColor','w','MarkerSize',3)
                end
            end
        end
        
        % add argo trajectory
        if argo
            for i=1:length(argos)
                dateb=datenum(argos(i).date)-datenum(dr);
                lonb=argos(i).lon;
                latb=argos(i).lat;
                mrub=find(dateb(dateb<=0)==max(dateb(dateb<=0)));
                ruban=find(dateb<=0 & dateb>-nba &...
                    diff([dateb;10000]) > 0.5 &...
                    lonb<maxlon & lonb>minlon & latb<maxlat & latb>minlat);
                if ~isempty(ruban)
                    L=[length(ruban(:)) i];
                    m_plot(lonb(ruban(:)),latb(ruban(:)),'-','color',[1 0 1],'linewidth',1)
                    m_plot(lonb(ruban(:)),latb(ruban(:)),'dk','markerfacecolor',[1 0 1],...
                        'markersize',3)
                    m_plot(lonb(ruban(end)),latb(ruban(end)),'dk','markerfacecolor',[.5 0 .5],...
                         'markersize',3)
                    %m_text(lonb(ruban(1)),latb(ruban(1))+0.1,argos(i).name)
                end
                if dateb(mrub)>-1
                    m_plot(lonb(mrub(end)),latb(mrub(end)),'dk','markerfacecolor',[.5 0 .5],...
                        'markersize',5)

                end
            end
            m_text(maxlon-5,maxlat+0.4,['ARGO last ',num2str(nba),'days (pink diamonds)'],'FontSize',6)
        end
        
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
            m_text(maxlon-5,maxlat+0.2,'EDDIES last 30days trajectories (grey lines)','FontSize',6)
            if strcmp(named,'BWE')
                m_text(5.65,44,'1 m/s','fontsize',6,'fontweight','bold')
            elseif strcmp(named,'LEV')
                m_text(29.4,37.25,'1 m/s','fontsize',6,'fontweight','bold')
            end
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
            set(hfig,'PaperPosition',[0,-1,10,6.5])
            set(hfig,'PaperSize',[8,6.5])
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
                mrub=find(dateb(dateb<=0)==max(dateb(dateb<=0)));
                ruban=find(dateb<=0 & dateb>-nba &...
                    diff([dateb;10000]) > 0.5 &...
                    lonb<maxlon & lonb>minlon & latb<maxlat & latb>minlat);
                if ~isempty(ruban)
                    m_plot(lonb(ruban(1)),latb(ruban(1)),'dk','markerfacecolor',[1 0 1],...
                        'markersize',6)
                    m_text(lonb(ruban(1)),latb(ruban(1))+0.15,num2str(argos(i).name),...
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
    
    crop([path_out,'/figures/',named{1},'/',num2str(n+1999)])
    crop([path_out,'/figures/',named{1},'_ID/',num2str(n+1999)])
    
    if movie, close(M); end
        
end

end

if export_out

%% export for Alex
% load tracking result
load([path_out,'eddy_tracks_',year])

INDA=[];stpA=[];LifetimeA=[];TminA=[];RmaxA=[];VmaxA=[];detamaxA=[];ellip_maxA=[];VortMA=[];
INDC=[];stpC=[];LifetimeC=[];TminC=[];RmaxC=[];VmaxC=[];detamaxC=[];ellip_maxC=[];VortMC=[];

% add eddy tracking (use get_aviso_nrt.sh in DATA/nrt -> mk_nc_along.m)
for i=1:length(tracks)
    dura = tracks(i).step(end)-tracks(i).step(1)+1;
    if tracks(i).type(1)==1
        INDC = [INDC;ones(length(tracks(i).step),1)*i];
        stpC = [stpC;tracks(i).step];
        LifetimeC = [LifetimeC;ones(length(tracks(i).step),1)*dura];
        TminC = [TminC;tracks(i).tau1];
        RmaxC = [RmaxC;tracks(i).rmax1];
        VmaxC = [VmaxC;tracks(i).velmax1];
        detamaxC = [detamaxC;tracks(i).deta1];
        ellip_maxC = [ellip_maxC;tracks(i).ellip1];
        VortMC = [VortMC;tracks(i).vortM1];
    else
        INDA = [INDA;ones(length(tracks(i).step),1)*i];
        stpA = [stpA;tracks(i).step];
        LifetimeA = [LifetimeA;ones(length(tracks(i).step),1)*dura];
        TminA = [TminA;tracks(i).tau1];
        RmaxA = [RmaxA;tracks(i).rmax1];
        VmaxA = [VmaxA;tracks(i).velmax1];
        detamaxA = [detamaxA;tracks(i).deta1];
        ellip_maxA = [ellip_maxA;tracks(i).ellip1];
        VortMA = [VortMA;tracks(i).vortM1];
    end
end

% plot distribution
load([path_out,'Atlas_main_features_AC'])
cut=15;
tps=0:15:1500;% in days
indA = LifetimeA>cut;
indC = LifetimeC>cut;
NLtA = histc(LifetimeA(indA),tps);
NLtC = histc(LifetimeC(indC),tps);
CSNA = cumsum(NLtA(end:-1:1));
CSNC = cumsum(NLtC(end:-1:1));

close all
hfig=figure('visible','on');
set(hfig,'Position',[0 0 1000 500])
set(gcf,'color','w'); % set figure background to white
subplot(1,2,1)
semilogy(tps,NLtC,'color',[.9 0 0],'linewidth',2)
hold on
semilogy(tps,NLtA,'color',[0 0 .9],'linewidth',2)
xlim([0 600])
xlabel('Life Time (days)','fontweight','bold','Fontsize',12)
ylabel('Number','fontweight','bold','Fontsize',12)
grid on
legend('Cyclones','Anticyclones')
subplot(1,2,2)
hold on
plot(tps,NLtC./NLtA,'color',[.7 0 .7],'linewidth',2)
plot(tps,CSNC(end:-1:1)./CSNA(end:-1:1),'color',[.5 .5 .5],'linewidth',2)
plot([0 1500],[1 1],'-k','linewidth',.5)
xlim([0 600])
xlabel('Life Time (days)','fontweight','bold','Fontsize',12)
ylabel('Cyclones / Anticyclones','fontweight','bold','Fontsize',12)
legend('15days bins','Cumulating longer eddies')
grid on
box on
%
set(hfig,'PaperPosition',[0,0,10,5])
print(hfig,[path_out,'figures/Distribution_Lifetime'],'-dpng','-r150')
crop([path_out,'/figures/'])

% Export Atlas linearised
% TA=table(INDA,stpA,LifetimeA,round(TminA),round(RmaxA*100)/100,round(VmaxA*100),...
%     round(detamaxA*100),round(ellip_maxA*100)/100,round(VortMA*10^7)/10^7,...
%     'VariableNames',{'Anticyclone_indice','step','Life_Time','Turnover_time_in_days',...
%         'Rmax_in_km','Vmax_in_cm_per_s','delta_ssh_in_cm','ellipticity','Maximal_vorticity'});
% writetable(TA,[path_out,'Atlas_main_features_A.dat'])
% %
% TC=table(INDC,stpC,LifetimeC,round(TminC),round(RmaxC*100)/100,round(VmaxC*100),...
%     round(detamaxC*100),round(ellip_maxC*100)/100,round(VortMC*10^7)/10^7,...
%     'VariableNames',{'Cyclone_indice','step','Life_Time','Turnover_time_in_days',...
%         'Rmax_in_km','Vmax_in_cm_per_s','delta_ssh_in_cm','ellipticity','Maximal_vorticity'});
% writetable(TA,[path_out,'Atlas_main_features_C.dat'])
%
save([path_out,'Atlas_main_features_AC'],'INDA','INDC','stpA','stpC',...
    'LifetimeA','LifetimeC','TminA','TminC','RmaxA','RmaxC','VmaxA','VmaxC',...
    'detamaxA','detamaxC','ellip_maxA','ellip_maxC','VortMA','VortMC')

end

