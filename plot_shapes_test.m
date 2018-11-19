% plot_tracking_test.m
%
% Make pdf with the shapes and centers
%
%=========================

start
clear; clc;

movie=0;
h_coast=1;
bath=1;

%----------------------------------------
% source of data driving the netcdf format
source = 'AVISO';

%----------------------------------------
% domaine
config = 'MED';
config = 'PRO2017_adt';

% Definition of the parameters specific for the experiment is needed 
run(['keys_sources_',source,'_',config])
load('param_eddy_tracking')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation vorticity, Okubo-Weiss and LNAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path_out2='/home/blevu/Resultats/AVISO/ALG/adt_2013/tests_v2/';
name='2000_2015';
name='2015';

%% Plot figrue

if movie || plo

% make figure's dir
% system(['mkdir ',path_out2,'figures']);
% system(['mkdir ',path_out2,'figures/pdf']);

% Load eddy detection result
%load([path_out,'fields_inter_2013'])
load([path_out,'eddy_centers_',name])
load([path_out,'eddy_shapes_',name])

% load bathymetry
if bath
    load('bathy/bathy_med')
end

% load high coastal def
if h_coast
    coast=csvread('bathy/new_bathy.csv');
end

% load grid
eval(['[x,y,mask,~,~,~] = load_fields_',source,'(1,resol,deg);'])

% set boundaries
minlat=min(y(:));
maxlat=max(y(:));
minlon=min(x(:));
maxlon=max(x(:));

% set boundaries
minlat=36;
maxlat=41;
minlon=-1;
maxlon=10;
% east
minlat=30;
maxlat=36;
minlon=20;
maxlon=32;
% IE04
minlat=31.5;
maxlat=34;
minlon=24;
maxlon=28;

% set the limit of eddies life
dayi=[2005 01 01 12 0 0];
dayi=[2014 01 01 12 0 0];
stepF=365;
st=[1 stepF];
dstp=1827;
dstp=365;

% set the limit of eddies life
cut=0;

m_proj('Mercator','lat',[minlat maxlat],'lon',[minlon maxlon]);

close all

for n=1:length(st)-1
    if movie
    % prepare fig
        hfig=figure('visible','off');
        set(hfig,'Position',[0 0 1000 600])
        M=VideoWriter([path_out2,'figures/shapes_',num2str(st(n)),'_',num2str(st(n+1)),'_new2.avi']);
        M.FrameRate = 1;
        open(M);
    end
    for day=1%st(n):st(n+1)
        if ~movie
            figure, close
            hfig=figure('visible','on');
            set(hfig,'Position',[0 0 1000 600])
        end
    % reference day at 12:00
        dr = datevec(datenum(dayi)+day-1);
    % pcolor with wwh
        CD01 = [centers0(day).x;centers0(day).y];
        CD02 = [centers0(day).x;centers0(day).y];
        CD11 = [centers(day).x;centers(day).y];
        CD12 = [centers(day).x;centers(day).y];
        CD21 = [centers2(day).x1;centers2(day).y1];
        CD22 = [centers2(day).x1;centers2(day).y1];
        eval(['[x,y,mask,u1,v1,ssh1] = load_fields_',source,'(day+dstp,1,deg);'])
        %eval(['[xi,yi,maski,~,~,~] = load_fields_',source,'(day+dstp,resol,deg);'])
        %L1=detection_fields(day).LNAM;
        %L2=L1;
        %LOW=detection_fields(day).LOW;
        %OW=detection_fields(day).OW;
        %vort=detection_fields(day).vort;
        %L1(OW>0)=nan;
        %L2(LOW>0)=nan;
        %maski(maski==0)=nan;
        xnan=x;
        ynan=y;
        vel=sqrt(u1.^2+v1.^2);
        xnan(isnan(vel) | vel<0.05) = nan;
        ynan(isnan(vel) | vel<0.05) = nan;
        %m_contourf(x,y,ssh1*100,50)
        %m_contourf(xi,yi,L2.*maski)
        shading flat
        colorbar
        caxis([-30 15])
        %caxis([-1 1])
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
        %ncol=50;
        %color=flipud(cbrewer('div','RdBu',ncol));
        %colormap(color(4:ncol-4,:,:))
        set(gcf,'color','w'); % set figure background to white
        % plot centers
        ind=find(CD01(1,:)<maxlon & CD01(1,:)>minlon & CD01(2,:)>minlat & CD01(2,:)<maxlat);
        m_plot(CD01(1,ind),CD01(2,ind),'xk','MarkerSize',4)
        ind=find(CD11(1,:)<maxlon & CD11(1,:)>minlon & CD11(2,:)>minlat & CD11(2,:)<maxlat);
        m_plot(CD11(1,ind),CD11(2,ind),'ok','MarkerFaceColor','r','MarkerSize',4)
        ind=find(CD21(1,:)<maxlon & CD21(1,:)>minlon & CD21(2,:)>minlat & CD21(2,:)<maxlat);
        m_plot(CD21(1,ind),CD21(2,ind),'ok','MarkerFaceColor','k','MarkerSize',4)
        for i=1:length(CD01)
            if CD01(1,i)<maxlon && CD01(1,i)>minlon && CD01(2,i)>minlat && CD01(2,i)<maxlat
                %m_text(CD01(1,i),CD01(2,i),['  ',num2str(i)],...
                 %   'color',[0 0 0],'FontSize',8,'fontWeight','bold')
            end
        end
        for i=1:length(CD11)
            if CD11(1,i)<maxlon && CD11(1,i)>minlon && CD11(2,i)>minlat && CD11(2,i)<maxlat
                %m_text(CD11(1,i),CD11(2,i),['  ',num2str(i)],...
                 %   'color',[0 0 0],'FontSize',8,'fontWeight','bold')
            end
        end
        for i=1:length(CD21)
            if CD21(1,i)<maxlon && CD21(1,i)>minlon && CD21(2,i)>minlat && CD21(2,i)<maxlat
                %m_text(CD21(1,i),CD21(2,i),['  ',num2str(i)],...
                   % 'color',[0 0 0],'FontSize',8,'fontWeight','bold')
            end
        end
        % plot shapes
        for i=1:length(shapes1(day).xy)
            if CD21(1,i)<maxlon && CD21(1,i)>minlon && CD21(2,i)>minlat && CD21(2,i)<maxlat
                lonlat1=shapes1(day).xy{i};
                lonlat2=shapes2(day).xy{i};
                lonlat3=shapes1(day).xy_end{i};
                type=centers2(day).type(i);
                calcul=warn_shapes2(day).calcul_curve(i);
                large1=warn_shapes2(day).large_curve1(i);
                if type==-1
                    col=[.1 .1 .9]; % anticyclone
                else
                    col=[.9 .1 .1]; % cyclone
                end
                if calcul==1
                    col(2)=.5; % calculated contour
                end
                if large1==0
                    wid='-'; % largest eddy
                else
                    wid='--'; % true eddy
                end
                if ~isempty(lonlat3)
                        m_plot(lonlat3(1,:),lonlat3(2,:),'--k','linewidth',1.5)
                end
                if ~isnan(lonlat1)
                    m_plot(lonlat1(1,:),lonlat1(2,:),wid,'color',col,'linewidth',2)
                end
                if ~isnan(lonlat2)
                    m_plot(lonlat2(1,:),lonlat2(2,:),'color',[.1 .9 .1],'linewidth',2) % double eddy
                end
            end
        end
        hold off
        % add legend on the map and information
        title(['AVISO from ',datestr(dayi,'dd mmm yyyy'),' for eddies living more than ',num2str(cut),' days'])
        %m_text(minlon+0.35,maxlat-0.3,[datestr(dr,1),' / degradation = ',num2str(deg)],...
        %    'FontSize',10,'fontWeight','bold')
        m_text(5.5,36.3,[datestr(dr,1),' / <gama> = ',num2str(0.8/deg)],...
            'FontSize',10,'fontWeight','bold')
        m_text(minlon+0.5,minlat+0.5,[datestr(dr,1)],...
            'FontSize',10,'fontWeight','bold')
        m_text(maxlon,maxlat+0.15,'ADT (cm)','FontSize',8,'fontWeight','bold')
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
            print(hfig,[path_out,'figures/pdf/tracking_',datestr(dr,1)],'-dpdf','-r0')
        end
    end

    if movie, close(M); end
    
end

end

