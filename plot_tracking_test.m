% plot_tracking.m
%
% Make a movie with tracks
%
%=========================

start
clear; clc;

movie=1;
h_coast=1;
bath=1;
tracking=1;
shaping=0;

%----------------------------------------
% source of data driving the netcdf format
source = 'AVISO';

%----------------------------------------
% domaine
dom = 'MED';

name='2013';

% Definition of the parameters specific for the experiment is needed 
%mod_eddy_params(['keys_sources_',source,'_',dom])
run(['keys_sources_',source,'_',dom])
load('param_eddy_tracking')

%% plot a movie or pdf
if movie || plo

% make figure's dir
system(['mkdir ',path_out,'figures']);
system(['mkdir ',path_out,'figures/pdf']);

% load tracking result
%load([path_out,'eddy_centers_',name])  
%load([path_out,'eddy_shapes_',name])  
load([path_out,'eddy_tracks_',name])  

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

% set area
minlat=min(y(:));
maxlat=max(y(:));
minlon=min(x(:));
maxlon=max(x(:));
minlat=36;
maxlat=41;
minlon=0;
maxlon=10;
%minlat=36;
%maxlat=41;
%minlon=5;
%maxlon=9.5;

% set period
dayi=[2013 01 01 12 0 0];
stepF=365;
st=[1 stepF];
dstp=4749;

% set the limit of eddies life
cut=0;

m_proj('Mercator','lat',[minlat maxlat],'lon',[minlon maxlon]);

close all

for n=1:length(st)-1
    if movie
    % prepare fig
        hfig=figure('visible','off');
        set(hfig,'Position',[0 0 1000 600])
        M=VideoWriter([path_out,'figures/test_1_',num2str(stepF),'_new3.avi']);
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
        eval(['[x,y,mask,u1,v1,ssh1] = load_fields_',source,'(day+dstp,1,deg);'])
        xnan=x;
        ynan=y;
        vel=sqrt(u1.^2+v1.^2);
        xnan(isnan(vel) | vel<0.05) = nan;
        ynan(isnan(vel) | vel<0.05) = nan;
        m_contourf(x,y,ssh1*100,50)
        %m_pcolor(lon,lat,LNAM.*mask)
        shading flat
        colorbar
        caxis([-20 25])
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
        tracks=tracks2;
        for i=1:length(tracks)
            ind=[];
            ind=find(tracks(i).step==day+dstp,1);
            if ~isempty(ind)
                CD = [(tracks(i).x1)';(tracks(i).y1)'];
                %CD = [(tracks(i).xbary1)';(tracks(i).ybary1)'];
                dura = tracks(i).step(ind)-tracks(i).step(1)+1;
                if tracks(i).step(end)-tracks(i).step(1)+1>=cut
                    m_plot(CD(1,1:ind),CD(2,1:ind),'-','color',[0.4 0.4 0.4],'linewidth',2)
                    m_plot(CD(1,1),CD(2,1),'sk','MarkerFaceColor','k','MarkerSize',4)
                    m_plot(CD(1,ind),CD(2,ind),'ok','MarkerFaceColor','k',...
                        'MarkerSize',round(dura/60+4))
                    %
                    lonlat3=tracks(i).shapes3{ind};
                    m_plot(lonlat3(1,:),lonlat3(2,:),'--k','linewidth',1.5)
                    %
                    lonlat1=tracks(i).shapes1{ind};
                    if ~isnan(lonlat1)
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
                    end
                    %
                    lonlat2=tracks(i).shapes2{ind};
                    if ~isnan(lonlat2)
                        m_plot(lonlat2(1,:),lonlat2(2,:),'-','color',[.1 .9 .1],'linewidth',2) % double eddy
                        list=unique(tracks2(i).interaction(1:ind));
                        list=list(~isnan(list) & list~=tracks2(i).interaction(ind))';
                        for j=list
                            ind11=find(tracks2(i).interaction==j,1,'first');
                            ind12=find(tracks2(i).interaction==j,1,'last');
                            ind21=find(tracks2(j).interaction==i,1,'first');
                            ind22=find(tracks2(j).interaction==i,1,'last');
                            CD2=[tracks2(i).x1(ind11),tracks2(i).x1(ind12);tracks2(i).y1(ind11),tracks2(i).y1(ind12)];
                            indj=find(tracks2(j).step<=day+dstp);
                            CD1 = [(tracks2(j).x1(indj))';(tracks2(j).y1(indj))'];
                            if ~isempty(ind21)
                            if nanmax(tracks2(j).split(ind21:ind22))==1 && nanmax(tracks2(j).merge(ind21:ind22))==1
                                CD1 = [CD2(:,1),CD1,CD2(:,end)];
                                col= [1 1 1]; % merging + splitting
                                per= [2 length(CD1)-1];
                            elseif nanmax(tracks2(j).split(ind21:ind22))==1
                                CD1 = [CD2(:,1),CD1];
                                col= [.6 0 .6]; % splitting
                                per= [2 length(CD1)];
                            elseif nanmax(tracks2(j).merge(ind21:ind22))==1
                                CD1 = [CD1,CD2(:,end)];
                                col=[0 .6 .6]; % merging
                                per=[1 length(CD1)-1];
                            end
                            if nanmax(tracks2(j).split(ind21:ind22))==1 || nanmax(tracks2(j).merge(ind21:ind22))==1
                                m_plot(CD1(1,:),CD1(2,:),'-','color',col,'linewidth',2)
                                m_plot(CD1(1,per),CD1(2,per),'dk','MarkerFaceColor','k','MarkerSize',4)
                            end
                            end
                        end
                        list=unique(tracks2(i).interaction2(1:ind));
                        list=list(~isnan(list) & list~=tracks2(i).interaction2(ind))';
                        for j=list
                            ind11=find(tracks2(i).interaction2==j,1,'first');
                            ind12=find(tracks2(i).interaction2==j,1,'last');
                            ind21=find(tracks2(j).interaction==i,1,'first');
                            ind22=find(tracks2(j).interaction==i,1,'last');
                            CD2=[tracks2(i).x1(ind11),tracks2(i).x1(ind12);tracks2(i).y1(ind11),tracks2(i).y1(ind12)];
                            indj=find(tracks2(j).step<=day+dstp);
                            CD1 = [(tracks2(j).x1(indj))';(tracks2(j).y1(indj))'];
                            if ~isempty(ind21)
                            if nanmax(tracks2(j).split(ind21:ind22))==1 && nanmax(tracks2(j).merge(ind21:ind22))==1
                                CD1 = [CD2(:,1),CD1,CD2(:,end)];
                                col= [1 1 1]; % merging + splitting
                                per= [2 length(CD1)-1];
                            elseif nanmax(tracks2(j).split(ind21:ind22))==1
                                CD1 = [CD2(:,1),CD1];
                                col= [.6 0 .6]; % splitting
                                per= [2 length(CD1)];
                            elseif nanmax(tracks2(j).merge(ind21:ind22))==1
                                CD1 = [CD1,CD2(:,end)];
                                col=[0 .6 .6]; % merging
                                per=[1 length(CD1)-1];
                            end
                            if nanmax(tracks2(j).split(ind21:ind22))==1 || nanmax(tracks2(j).merge(ind21:ind22))==1
                                m_plot(CD1(1,:),CD1(2,:),'-','color',col,'linewidth',2)
                                m_plot(CD1(1,per),CD1(2,per),'dk','MarkerFaceColor','k','MarkerSize',4)
                            end
                            end
                        end
                    end                    
                    %
                    if CD(1,ind)<maxlon && CD(1,ind)>minlon && CD(2,ind)>minlat && CD(2,ind)<maxlat
                        if tracks(i).split(ind)==1
                            m_text(CD(1,ind),CD(2,ind)+0.1,'split',...
                                'color',[0 0 0],'FontSize',8,'fontWeight','bold')
                        end
                        if tracks(i).merge(ind)==1
                            m_text(CD(1,ind),CD(2,ind)-0.1,'merge',...
                                'color',[0 0 0],'FontSize',8,'fontWeight','bold')
                        end
                        %m_text(CD(1,ind),CD(2,ind),['  ',num2str(dura)],...
                          %  'color',[0 0 0],'FontSize',round(dura/60+4),'fontWeight','bold')
                        m_text(CD(1,ind),CD(2,ind)-.2,[num2str(i)],...
                            'color',[.5 .5 .5],'FontSize',8,'fontWeight','bold')
                    end
                else
                    m_plot(CD(1,ind),CD(2,ind),'ok',...
                        'MarkerFaceColor','w','MarkerSize',4)
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
        m_text(9.6,41.2,'ADT (cm)','FontSize',8,'fontWeight','bold')
        %m_text(9.5,40.7,'ADT (cm)','FontSize',8,'fontWeight','bold')
         %prepare the print in a pdf
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

