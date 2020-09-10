% plot contours and center on a specific day
start
clear; clc;
source = 'AVISO';
config = 'MED_adt';
day = 164; % time step
dayi=[2019 05 01 12 0 0];% from data file
dr = datevec(datenum(dayi)+day-1);

% load parameters, centers and shapes results
run(['keys_sources_',source,'_',config])
load('param_eddy_tracking'); load([path_out,'eddy_centers']); load([path_out,'eddy_shapes'])
CD21 = [centers2(day).x1;centers2(day).y1];

% load grid and fields at day
eval(['[x,y,mask,ssu,ssv,ssh] = load_fields_',source,'(day,1,deg);'])

% set boundaries and plot
minlat=36.5; maxlat=40.5; minlon=3.5; maxlon=9;
m_proj('Mercator','lat',[minlat maxlat],'lon',[minlon maxlon]);

% plot the day
close all
hfig=figure('visible','off');
set(hfig,'Position',[0 0 800 800])

% plot quiver velocities
xnan=x; ynan=y; vel=sqrt(ssu.^2+ssv.^2);
xnan(isnan(vel) | vel<0.05) = nan; ynan(isnan(vel) | vel<0.05) = nan;
m_quiver(xnan,ynan,ssu,ssv,2,'color',[0.3 0.3 0.3])

% grid and fancy the map
m_coast('patch',[.9 .9 .9]); m_grid('tickdir','in','linewidth',1,'linestyle','none');
set(gcf,'color','w'); % set figure background to white

% plot shapes
hold on
for i=1:length(shapes1(day).xy)
  if CD21(1,i)<maxlon && CD21(1,i)>minlon && CD21(2,i)>minlat && CD21(2,i)<maxlat
     lonlat1=shapes1(day).xy{i}; lonlat2=shapes2(day).xy{i}; lonlat3=shapes1(day).xy_end{i};
     type=centers2(day).type(i);
     if type==-1, 
       col=[.1 .1 .9]; % anticyclone
     else
       col=[.9 .1 .1]; % cyclone
     end
     if ~isempty(lonlat3)
       hl(3) = m_plot(lonlat3(1,:),lonlat3(2,:),'--k','linewidth',1.5); % last contour
     end
     if ~isnan(lonlat1)
              if type==-1
	  hl(1) = m_plot(lonlat1(1,:),lonlat1(2,:),'color',col,'linewidth',2); % characteristic contour
	else
	  hl(2) = m_plot(lonlat1(1,:),lonlat1(2,:),'color',col,'linewidth',2); % characteristic contour
	end
     end

     if ~isnan(lonlat2)
       hl(4) = m_plot(lonlat2(1,:),lonlat2(2,:),'color',[.1 .9 .1],'linewidth',2); % shared contour
     end
   end
end

% plot centers
ind=find(CD21(1,:)<maxlon & CD21(1,:)>minlon & CD21(2,:)>minlat & CD21(2,:)<maxlat);
m_plot(CD21(1,ind),CD21(2,ind),'ok','MarkerFaceColor','k','MarkerSize',4)
for i=1:length(CD21)
  if CD21(1,i)<maxlon && CD21(1,i)>minlon && CD21(2,i)>minlat && CD21(2,i)<maxlat
    m_text(CD21(1,i),CD21(2,i),['  ',num2str(i)],'color',[0 0 0],'FontSize',10,'fontWeight','bold')
  end
end
hold off

% add legend on the map and information and print out
title('AMEDA on geostrophic velocities from AVISO','Fontsize',14)
m_text(minlon+0.2,minlat+0.2,[datestr(dr,1)],'FontSize',10,'fontWeight','bold')
legend(hl,{'characteristic contour (AE)','characteristic contour (CE)','last contour','shared contour'})
set(hfig,'PaperPosition',[0,0,8,8]), set(hfig,'PaperSize',[8,8])
print(hfig,[path_out,'shapes_',datestr(dr,'yyyymmdd')],'-dpdf','-r100')

