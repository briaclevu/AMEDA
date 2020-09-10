% define configuration and the time
start, clear; clc;
source = 'AVISO'; config = 'MED_adt';
day = 164; % time step
dayi=[2019 05 01 12 0 0];% from data file
dr = datevec(datenum(dayi)+day-1);
cut=0; % filter the results or not (=0)
dstp=0;

% load parameters and tracks results
run(['keys_sources_',source,'_',config])
load('param_eddy_tracking'); load([path_out,'eddy_tracks2'])

% set boundaries and plot
minlat=36.5; maxlat=40.5; minlon=3.5; maxlon=9;
m_proj('Mercator','lat',[minlat maxlat],'lon',[minlon maxlon]);

% define a plot window and the map
close all, clear hl
hfig=figure('visible','off');
set(hfig,'Position',[0 0 800 800])

% load grid and fields at the step day
eval(['[x,y,mask,ssu,ssv,~] = load_fields_',source,'(day,1,deg);'])

% plot quiver velocities
xnan=x; ynan=y; vel=sqrt(ssu.^2+ssv.^2);
xnan(isnan(vel) | vel<0.05) = nan; ynan(isnan(vel) | vel<0.05) = nan;
m_quiver(xnan,ynan,ssu,ssv,2,'color',[0.3 0.3 0.3])

% grid and fancy the map
m_coast('patch',[.9 .9 .9]); m_grid('tickdir','in','linewidth',1,'linestyle','none');
set(gcf,'color','w'); % set figure background to white

hold on
% add eddy tracking
tracks=tracks2;
for i=1:length(tracks)
    ind=[];
    ind=find(tracks(i).step==day,1);
    if ~isempty(ind)
        CD = [(tracks(i).xbary1)';(tracks(i).ybary1)'];
        dura = tracks(i).step(ind)-tracks(i).step(1)+1;
        if tracks(i).step(end)-tracks(i).step(1)+1>=cut
            hl(5) = m_plot(CD(1,1:ind),CD(2,1:ind),'-','color',[0.4 0.4 0.4],'linewidth',2);
            m_plot(CD(1,1),CD(2,1),'sk','MarkerFaceColor','k','MarkerSize',4)
            m_plot(CD(1,ind),CD(2,ind),'ok','MarkerFaceColor','k',...
                'MarkerSize',round(dura/60+4))
            % plot last contour
            lonlat3=tracks(i).shapes3{ind};
            hl(3) = m_plot(lonlat3(1,:),lonlat3(2,:),'--k','linewidth',1.5);
            % plot characteristic contour
            lonlat1=tracks(i).shapes1{ind};
            if ~isnan(lonlat1)
                if tracks(i).type(1)==-1
                    col=[.1 .1 .9]; % anticyclone
                    hl(1) = m_plot(lonlat1(1,:),lonlat1(2,:),'color',col,'linewidth',2);
                else
                    col=[.9 .1 .1]; % cyclone
                    hl(2) = m_plot(lonlat1(1,:),lonlat1(2,:),'color',col,'linewidth',2);
                end
            end
            % plot shared contour
            lonlat2=tracks(i).shapes2{ind};
            if ~isnan(lonlat2)
                hl(4) = m_plot(lonlat2(1,:),lonlat2(2,:),'-','color',[.1 .9 .1],'linewidth',2); % double eddy
            end
            %
            if CD(1,ind)<maxlon && CD(1,ind)>minlon && CD(2,ind)>minlat && CD(2,ind)<maxlat
                if tracks(i).split(ind)==1
                    m_text(CD(1,ind)+0.2,CD(2,ind),'split',...
                        'color',[0 0 0],'FontSize',8,'fontWeight','bold')
                end
                if tracks(i).merge(ind)==1
                    m_text(CD(1,ind)+0.2,CD(2,ind),'merge',...
                        'color',[0 0 0],'FontSize',8,'fontWeight','bold')
                end
                %m_text(CD(1,ind),CD(2,ind),['  ',num2str(dura)],...
                  %  'color',[0 0 0],'FontSize',round(dura/60+4),'fontWeight','bold')
                m_text(CD(1,ind),CD(2,ind)+.1,[num2str(i)],...
                    'color',[0 0 0],'FontSize',10,'fontWeight','bold')
            end
        end
    end
end

% add legend on the map and information
title('AMEDA tracks on geostrophic velocities from AVISO','Fontsize',14)
m_text(minlon+0.2,minlat+0.2,[datestr(dr,1)],'FontSize',10,'fontWeight','bold')
legend(hl,{'Anticyclone','Cyclone','last contour','Interaction contour','Barycenter tracks'})
set(hfig,'PaperPosition',[0,0,8,8]), set(hfig,'PaperSize',[8,8])
print(hfig,[path_out,'tracking_',datestr(dr,'yyyymmdd')],'-dpdf','-r100')

