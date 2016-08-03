% MAIN_eddy_tracking.m
%
%   MAIN_eddy_tracking is the main function of the eddy detection and
%   tracking package. It returns position of the centers, dimensions and 
%   tracks of the eddies detected from the time series of a 2-D velocity 
%   field.
%   It gives also an history of the splitting and merging events.
%
%   The algortihm subroutines:
%
%   - param_eddy_tracking sets user defined paths and parameters:
%     nc_u nc_v nc_dim b bx r path_in path_out periodic criteres
%     Users should modify param_eddy_tracking.m according to their 
%     settings.
%     IMPORTANT: See also the documentation from param_eddy_tracking.m for 
%     a description of the format requirements for the input files.
%
%   - mod_eddy_centers returns a structure array with the position of the
%     detected eddy centers.
%
%   - mod_eddy_shapes computes features for the detected eddy centers.
%
%   - mod_eddy_tracks computes eddy tracks using the detected centers.
%
%   The package was developed using Matlab 8 R2014a, and it requires few 
%   functions from the statistics toolbox (nanmean,...) and from Buehren
%   routine for assigment (assigmentoptimal).
%   It uses also netcdf high-level function (ncread).

%-------------------------
%
%   The routines are originally from the Nencioli et al.'s method.
%
%   Nencioli F., Dong C., Dickey T., Washburn L. and McWilliams J.,
%   "A Vector Geometry Based Eddy Detection Algorithm and Its Application 
%    to a High-resolution Numerical Model Product and High-frequency Radar 
%    Surface Velocities in the Southern California Bight".
%    2010, JAOT, Vol. 27, No. 3, pp. 564-579.
%
%   The original version has been modify using for the first time the
%   Local Normalized Angular Momentum (LNAM) added to the Nencioli et al.'s
%   method.
%
%   Nadia Mkhinini, Andre Luis Santi Coimbra, Alexandre Stegner,
%   Thomas Arsouze, Isabelle Taupier-Letage and Karine Béranger.
%   "Long-lived mesoscale eddies in the eastern Mediterranean Sea:
%    Analysis of 20 years of AVISO geostrophic velocities".
%   Journal of Geophysical Research. Oceans, Wiley-Blackwell, 2014,
%   pp.10.1002/2014JC010176.
%
%   The present method keep the structure and names of the original routine
%   but uses only the LNAM and Okubo-Weiss criteria for the centers
%   detection. It is the purpose of a paper in preparation.
%
%-------------------------
%   Ver. 3.2 2015 by Le Vu et al.
%   Ver. 3.1 2014 by Mkhinini et al. and LMD from from Nencioli et al routines's
%-------------------------
%
%=========================

clear; clc;

% Definition of the parameters specific for the experiment is needed 
deg=1; param_eddy_tracking_AVISO

global runname
global path_out

proc=0;
plo=0; % figure if ==1
movie=1;
h_coast=1;
bath=1;

%% Compute the entire serie
upd = 0;

if proc

% 1: get fields and grids u, v and ssh and interpolte them by a factor 'res'
load_fields_AVISO(nc_dim,nc_u,nc_v,nc_ssh,b,bx,res,deg);

% 2: Creation of interpolated fields used for detection used for each day
mod_fields_inter(0)

if res==1
    % non interpolated fields is interpoleted fields
    copyfile([path_out,'fields_inter',runname,'.mat'],[path_out,'fields',runname,'.mat'])
else
    % Creation of non interpolated fields used for detection used for each day
    mod_fields(0)
end

% 3: Detection of eddy centers
disp(['Find center for K = ',num2str(K)])
mod_eddy_centers(K,0)

% 4: Definition of eddy contours
disp('Find shapes')
mod_eddy_shapes(0)

% 5: Tracking of the eddies
cut_off=0;% save minimal duration (0=tau(n))
Dt=10;% searching delay tolerance
disp(['Tracks eddies at cut_off ',num2str(cut_off),' allowing a delay of ',num2str(Dt),' days']) 
mod_eddy_tracks(cut_off,Dt,T,dps,0)

end

%% plot a movie or pdf
if movie || plo

% make figure's dir
system(['mkdir ',path_out,'../figures']);
system(['mkdir ',path_out,'../figures/pdf']);

% load data fields
load([path_out,'fields',name,'.mat'])

% Load eddy detection result
load([path_out,'eddy_centers',name])
load([path_out,'eddy_shapes',name])

% load tracking result
load([path_out,'eddy_tracks',name])  

% load bathymetry
if bath
    load('bathy/bathy_med')
end

% load high coastal def
if h_coast
    coast=csvread('bathy/new_bathy.csv');
end

% set area
minlat=min(y(:));
maxlat=max(y(:));
minlon=min(x(:));
maxlon=max(x(:));

% set periodi
st=[1 size(u,3)];
dayi=[2013 01 01 12 0 0];

% set the limit of eddies life
cut_off=1;
cut=0;

m_proj('Mercator','lat',[minlat maxlat],'lon',[minlon maxlon]);

close all

for n=1%:length(st)-1
    if movie
    % prepare fig
        hfig=figure('visible','off');
        set(hfig,'Position',[0 0 1000 600])
        M=VideoWriter([path_out,'../figures/',sshname,num2str(st(n)),'_',num2str(st(n+1)),'_sm2.avi']);
        M.FrameRate = 1;
        open(M);
    end
    for day=st(n):st(n+1)
        if ~movie
            figure, close
            hfig=figure('visible','on');
            set(hfig,'Position',[0 0 1000 600])
        end
    % reference day at 12:00
        dr = datevec(datenum(dayi)+day-1);
    % pcolor with wwh
        ssh1=squeeze(ssh(:,:,day));
        u1=squeeze(u(:,:,day));
        v1=squeeze(v(:,:,day));
        xnan=x;
        ynan=y;
        vel=sqrt(u1.^2+v1.^2);
        xnan(isnan(vel) | vel<0.05) = nan;
        ynan(isnan(vel) | vel<0.05) = nan;
        m_contourf(x,y,ssh1*100,50)
        %m_pcolor(lon,lat,LNAM.*mask)
        shading flat
        colorbar
        if strcmp(sshname,'adt_')
            caxis([-20 25])
        else
            caxis([-15 10])
        end
        % add quiver velocities
        hold on
        m_quiver(xnan,ynan,u1,v1,3,'color',[0.3 0.3 0.3])
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
            m_coast('color','k');
        end
        % grid and fancy the map
        m_grid('tickdir','in','linewidth',1,'linestyle','none');
        colormap(soft(40,0.3))
        %colormap(cbrewer('seq','Oranges',20))
        set(gcf,'color','w'); % set figure background to white
        
        % add eddy tracking (use get_aviso_nrt.sh in DATA/nrt -> mk_nc_along.m)
        for i=1:length(tracks)
            ind=[];
            ind=find(tracks(i).step==day,1);
            if ~isempty(ind)
                CD = [(tracks(i).x1)';(tracks(i).y1)'];
                %CD = [(tracks(i).xbary1)';(tracks(i).ybary1)'];
                dura = tracks(i).step(ind)-tracks(i).step(1)+1;
                if tracks(i).step(end)-tracks(i).step(1)+1>=cut
                    m_plot(CD(1,1:ind),CD(2,1:ind),'-','color',[0.4 0.4 0.4],'linewidth',2)
                    if ~isnan(tracks(i).shapes2{ind})
                        lonlat2=tracks(i).shapes2{ind};
                        m_plot(lonlat2(1,:),lonlat2(2,:),'color',[.1 .9 .1],'linewidth',2) % double eddy
                    end
                    lonlat1=tracks(i).shapes1{ind};
                    if tracks(i).type(1)==-1
                        col=[.1 .1 .9]; % anticyclone
                    else
                        col=[.9 .1 .1]; % cyclone
                    end
                    if tracks(i).calcul(ind)==1
                        col(2)=.5; % calculated contour
                    end
                    if tracks(i).large1(ind)==0
                        wid='-'; % largest eddy
                    else
                        wid='--'; % true eddy
                    end
                    m_plot(lonlat1(1,:),lonlat1(2,:),wid,'color',col,'linewidth',2)
                    m_plot(CD(1,1),CD(2,1),'sk','MarkerFaceColor','k','MarkerSize',4)
                    m_plot(CD(1,ind),CD(2,ind),'ok','MarkerFaceColor','k',...
                        'MarkerSize',round(dura/60+4))
                    if CD(1,ind)<maxlon && CD(1,ind)>minlon && CD(2,ind)>minlat && CD(2,ind)<maxlat
                        if tracks(i).split(ind)==1
                            m_text(CD(1,ind),CD(2,ind)+0.1,'split',...
                                'color',[0 0 0],'FontSize',8,'fontWeight','bold')
                        end
                        if tracks(i).merge(ind)==1
                            m_text(CD(1,ind),CD(2,ind)-0.1,'merge',...
                                'color',[0 0 0],'FontSize',8,'fontWeight','bold')
                        end
                        m_text(CD(1,ind),CD(2,ind),['  ',num2str(dura)],...
                            'color',[0 0 0],'FontSize',round(dura/60+4),'fontWeight','bold')
                        %m_text(CD(1,ind),CD(2,ind),['   ',num2str(i)],...
                         %   'color',[0 0 0],'FontSize',8,'fontWeight','bold')
                    end
                else
                    m_plot(CD(1,ind),CD(2,ind),'ok',...
                        'MarkerFaceColor','w','MarkerSize',4)
                end
            end
        end
        
        % merging/splitting
        for i=1:length(tracks)
            ind=[];
            ind=find(tracks(i).step==day,1);
            if ~isempty(ind) && tracks(i).step(end)-tracks(i).step(1)+1>=cut
                CD1 = [(tracks(i).x1)';(tracks(i).y1)'];
                list=unique(tracks(i).interaction(1:ind));
                list=list(~isnan(list) & list~=tracks(i).interaction(ind))';
                for j=list
                    ind11=find(tracks(i).interaction==j,1,'first');
                    ind12=find(tracks(i).interaction==j,1,'last');
                    ind21=find(tracks(j).interaction==i,1,'first');
                    ind22=find(tracks(j).interaction==i,1,'last');
                    CD2=[tracks(i).x1(ind11),tracks(i).x1(ind12);tracks(i).y1(ind11),tracks(i).y1(ind12)];
                    indj=find(tracks(j).step<=day);
                    CD = [(tracks(j).x1(indj))';(tracks(j).y1(indj))'];
                    if nanmax(tracks(j).split(ind21:ind22))==1 && nanmax(tracks(j).merge(ind21:ind22))==1
                        CD = [CD2(:,1),CD,CD2(:,end)];
                        col= [1 1 1]; % merging + splitting
                        per= [2 length(CD)-1];
                    elseif nanmax(tracks(j).split(ind21:ind22))==1
                        CD = [CD2(:,1),CD];
                        col= [.6 0 .6]; % splitting
                        per= [2 length(CD)];
                    elseif nanmax(tracks(j).merge(ind21:ind22))==1
                        CD = [CD,CD2(:,end)];
                        col=[0 .6 .6]; % merging
                        per=[1 length(CD)-1];
                    end
                    if nanmax(tracks(j).split(ind21:ind22))==1 || nanmax(tracks(j).merge(ind21:ind22))==1
                        m_plot(CD(1,:),CD(2,:),'-','color',col,'linewidth',2)
                        m_plot(CD(1,per),CD(2,per),'dk','MarkerFaceColor','k','MarkerSize',4)
                    end
                end
                list=unique(tracks(i).interaction2(1:ind));
                list=list(~isnan(list) & list~=tracks(i).interaction2(ind))';
                for j=list
                    ind11=find(tracks(i).interaction2==j,1,'first');
                    ind12=find(tracks(i).interaction2==j,1,'last');
                    ind21=find(tracks(j).interaction==i,1,'first');
                    ind22=find(tracks(j).interaction==i,1,'last');
                    CD2=[tracks(i).x1(ind11),tracks(i).x1(ind12);tracks(i).y1(ind11),tracks(i).y1(ind12)];
                    indj=find(tracks(j).step<=day);
                    CD = [(tracks(j).x1(indj))';(tracks(j).y1(indj))'];
                    if nanmax(tracks(j).split(ind21:ind22))==1 && nanmax(tracks(j).merge(ind21:ind22))==1
                        CD = [CD2(:,1),CD,CD2(:,end)];
                        col= [1 1 1]; % merging + splitting
                        per= [2 length(CD)-1];
                    elseif nanmax(tracks(j).split(ind21:ind22))==1
                        CD = [CD2(:,1),CD];
                        col= [.6 0 .6]; % splitting
                        per= [2 length(CD)];
                    elseif nanmax(tracks(j).merge(ind21:ind22))==1
                        CD = [CD,CD2(:,end)];
                        col=[0 .6 .6]; % merging
                        per=[1 length(CD)-1];
                    end
                    if nanmax(tracks(j).split(ind21:ind22))==1 || nanmax(tracks(j).merge(ind21:ind22))==1
                        m_plot(CD(1,:),CD(2,:),'-','color',col,'linewidth',2)
                        m_plot(CD(1,per),CD(2,per),'dk','MarkerFaceColor','k','MarkerSize',4)
                    end
                end
            end
        end
        hold off
        % add legend on the map and information
        if cut_off==0
            title(['AVISO from ',datestr(dayi,'dd mmm yyyy'),' for eddies living more than turnover time'])
        else
            title(['AVISO from ',datestr(dayi,'dd mmm yyyy'),' for eddies living more than ',num2str(cut_off),' days'])
        end
        %m_text(minlon+0.35,maxlat-0.3,[datestr(dr,1),' / degradation = ',num2str(deg)],...
        %    'FontSize',10,'fontWeight','bold')
        m_text(5.5,36.3,[datestr(dr,1)],...%,' / <gama> = ',num2str(0.8/deg)],...
            'FontSize',10,'fontWeight','bold')
        m_text(9.5,40.7,'ADT (cm)','FontSize',8,'fontWeight','bold')
        % prepare the print in a pdf
        if movie
            figinfo = hardcopy(hfig,'-dzbuffer','-r0');
            writeVideo(M,im2frame(figinfo));
            %frame = getframe(gcf);
            %writeVideo(M,frame);
        else
            set(hfig,'PaperPosition',[-1.3,-0.4,10,7])
            set(hfig,'PaperSize',[7.5,6.5])
            %set(hfig,'PaperPosition',[-1,-0.5,9,8])
            %set(hfig,'PaperSize',[7,7.5])
            print(hfig,[path_out,'../figure/pdf/tracking_',datestr(dr,1)],'-dpdf','-r0')
        end
    end

    if movie, close(M); end
    
end

end