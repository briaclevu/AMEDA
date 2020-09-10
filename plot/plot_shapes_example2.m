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
load('param_eddy_tracking'); load([path_out,'eddy_shapes'])

% following the script of the example 1
% plot the V-R profil
clear hl, close all
hfig=figure('visible','off');
set(hfig,'Position',[0 0 600 400])
hold on

% double eddy center 78 with 65
i=78; j=65;
hl(1) = plot(profil2(day).rmoy{i},profil2(day).vel{i},'.');
hl(2) = plot(profil2(day).rmoy{j},profil2(day).vel{j},'.k');
hl(3) = plot(shapes1(day).rmax(i),shapes1(day).velmax(i),'or','markerfacecolor','r');
plot(shapes1(day).rmax(j),shapes1(day).velmax(j),'or','markerfacecolor','r');
hl(4) = plot(shapes2(day).rmax(j),shapes2(day).velmax(j),'og','markerfacecolor','g');

% add legend on the map and information and print out
axis([0 100 0 .4])
ylabel('<V> (m.s^-^1)','Fontsize',12)
xlabel('<R> (km)','Fontsize',12)
title(['<V>-<R> profiles from streamlines around eddies ',num2str(i),' and ',num2str(j)],'Fontsize',14)
legend(hl,{['Profile from eddy ',num2str(i)],…
                 ['Profile from eddy ',num2str(j)],…
                 ['V_m_a_x and R_m_a_x for main centers'],…
                 'V_m_a_x and R_m_a_x for shared contour'},'Location','Southeast')
grid on
box on

% print out
set(hfig,'PaperPosition',[0,0,6,4]), set(hfig,'PaperSize',[6,4])
print(hfig,[path_out,'profile_',num2str(i),'_',num2str(j),'_',datestr(dr,'yyyymmdd')],'-dpdf','-r100')

