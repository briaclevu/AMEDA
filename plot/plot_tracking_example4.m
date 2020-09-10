% define configuration
start, clear; clc;
source = 'AVISO'; config = 'MED_adt';
dayi=[2019 05 01 12 0 0];% from data file
janp=246+datenum(dayi)-1;

% load parameters and tracks results
run(['keys_sources_',source,'_',config])
load('param_eddy_tracking'); load([path_out,'eddy_tracks2']), tracks=tracks2;

% longest eddy
n=76; m=645; % track indices
stepn = tracks(n).step+datenum(dayi)-1;
stepm = tracks(m).step+datenum(dayi)-1;

% plot features eddy tracks n
close all
hfig=figure('visible','off');
set(hfig,'Position',[0 0 1000 800])
set(gcf,'color','w'); % set figure background to white

% time series
% Rmax
subplot(2,2,1)
plot(stepn,tracks(n).rmax1,'color',[0 .2 .8],'linewidth',2)
ylim([0 80])
hold on, plot([janp janp],[0 80],'k','linewidth',1.5)
text(janp-5,5,'2019','HorizontalAlignment', 'right','Fontsize',12), text(janp+5,5,'2020','Fontsize',12)
datetick('x','dd/mm','keepticks')
box on, grid on
xlabel('date','FontSize',14), ylabel('R_m_a_x (km)','FontSize',12)
title(['Size Anticyclone Eddy ',num2str(n)],'FontSize',14)
subplot(2,2,2)
plot(stepm,tracks(m).rmax1,'color',[.8 .2 0],'linewidth',2)
ylim([0 80])
text(stepm(end)-5,5,'2019','HorizontalAlignment', 'right','Fontsize',12)
datetick('x','dd/mm','keepticks')
box on, grid on
xlabel('date','FontSize',14), ylabel('R_m_a_x (km)','FontSize',12)
title(['Size Cyclone Eddy ',num2str(m)],'FontSize',14)

% Velmax
subplot(2,2,3)
plot(stepn,tracks(n).velmax1,'color',[0 .2 .8],'linewidth',2)
ylim([0 .6])
hold on, plot([janp janp],[0 .6],'k','linewidth',1.5)
text(janp-5,.05,'2019','HorizontalAlignment', 'right','Fontsize',12), text(janp+5,.05,'2020','Fontsize',12)
datetick('x','dd/mm','keepticks')
box on, grid on
xlabel('date','FontSize',14), ylabel('V_m_a_x (m.s^-^1)','FontSize',12)
title(['Intensity Anticyclone Eddy ',num2str(n)],'FontSize',14)
subplot(2,2,4)
plot(stepm,tracks(m).velmax1,'color',[.8 .2 0],'linewidth',2)
ylim([0 .6])
text(stepm(end)-5,.05,'2019','HorizontalAlignment', 'right','Fontsize',12)
datetick('x','dd/mm','keepticks')
box on, grid on
xlabel('date','FontSize',14), ylabel('V_m_a_x (m.s^-^1)','FontSize',12)
title(['Intensity Cyclone Eddy ',num2str(m)],'FontSize',14)

% print figure
set(hfig,'PaperPosition',[0,0,10,8])
set(hfig,'PaperSize',[10,8])
print(hfig,[path_out,'serie_eddy_',num2str(n),'_',num2str(m)],'-dpdf','-r100')



