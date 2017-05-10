close all

if exist('fig_num','var') == 0
    fig_num = 1;
end

screen_size = get(0, 'ScreenSize');

GlacialGrey = [0.4 0.7 1];
GargantuanGreen = [0.1 0.5 0.2];
OrdinaryOrange = [0.8 0.4 0];
ViolentViolet = [0.6 0 0.7];
MagicMaroon = [0.65 0.32 0.35];
EcstaticEmerald = [0.17 0.52 0.5];
GorgeousGold = [1 0.85 0];
BohemeBlue = [0.1 0.1 0.45];
CandidCoral = [1 0.44 0.32];
GuardianGreen = [0 0.8 0];
LightGrey = [0.65 0.65 0.65];
DarkGrey = [0.4 0.4 0.4];
matlabblue = [0 0.447 0.741];
matlabred = [0.85 0.325 0.098];
matlabyellow = [0.929 0.694 0.125];

% effective emissivity of the sky
ema_eff = nan(size(T_air_ac));
for t=1:length(time)
ema_eff(t) = LW_in_ac_val(t)/(5.67*10^(-8)*T_air_ac(t)^4);
end

load LWsub_HM.mat

% relative model error for sub-canopy LW radiation
LW_in_bc_CLM_RelErr = nan(size(LW_in_bc_CLM));
LW_in_bc_SP_RelErr = nan(size(LW_in_bc_SP));
for t=1:length(time)
    LW_in_bc_CLM_RelErr(t) = (LW_in_bc_CLM(t) - LW_in_bc_val(t))/LW_in_bc_val(t) * 100;
    LW_in_bc_SP_RelErr(t) = (LW_in_bc_SP(t) - LW_in_bc_val(t))/LW_in_bc_val(t) * 100;
end

LWenh_CLM = nan(size(LW_in_bc_CLM));
LWenh_CLM_HM = nan(size(LW_in_bc_CLM_HM));
LWenh_SP = nan(size(LW_in_bc_SP));
LWenh_val = nan(size(LW_in_bc_val));
for t=1:length(time)
    LWenh_CLM(t) = LW_in_bc_CLM(t)/LW_in_ac_val(t);
    LWenh_CLM_HM(t) = LW_in_bc_CLM_HM(t)/LW_in_ac_val(t);
    LWenh_SP(t) = LW_in_bc_SP(t)/LW_in_ac_val(t);
    LWenh_val(t) = LW_in_bc_val(t)/LW_in_ac_val(t);
end

%-------------------------------  Figures  --------------------------------
% Meteorology
x = time(4225:end);
y1 = LW_in_ac_val(4225:end);
y2 = T_air_ac(4225:end)-273.15;
y3 = SW_in_ac_val(4225:end);
y4 = nan(size(time(4225:end)));
y5 = nan(size(time(4225:end)));
y6 = ema_eff(4225:end);
timeticks = [time(4009) time(4418) time(4765) time(5114) time(5461) time(5858) time(6205) time(6578) time(6925) datenum(2004,5,31,23,59,59)];
datelabel = {'15 Jan' '01 Feb' '15 Feb' '01 Mar' '15 Mar' '01 Apr' '15 Apr' '01 May' '15 May' '01 Jun'};

fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
[ax,h1,h2] = plotyy(x,y1,x,y2);
ax(3) = axes('yaxislocation','left','Color','none','XColor','k','YColor',matlabyellow,'box','off');
h3 = line(x,y3,'Parent',ax(3),'Color',matlabyellow);
ax(4) = axes('yaxislocation','right','Color','none','XColor','k','YColor',GargantuanGreen,'box','off');
h4 = line(x,y6,'Parent',ax(4),'Color',GargantuanGreen);
set(ax(1),'XTick',timeticks,'XTickLabel',datelabel);%set(ax(1),'XTick',[4418 5114 5858 6578]); datetick(ax(1),'x','dd mmm')
set(ax(2),'XTick',timeticks,'XTickLabel',datelabel);%set(ax(2),'XTick',[4418 5114 5858 6578]); datetick(ax(2),'x','dd mmm')
set(ax(3),'XTick',timeticks,'XTickLabel',datelabel);%set(ax(3),'XTick',[4418 5114 5858 6578]); datetick(ax(3),'x','dd mmm')
set(ax(4),'XTick',timeticks,'XTickLabel',datelabel);%set(ax(4),'XTick',[4418 5114 5858 6578]); datetick(ax(4),'x','dd mmm')
set(ax(1),'Box','off')
set(ax(2),'Box','off')
set(ax(3),'Box','off')
set(ax(4),'Box','off')
xlabel('Alptal observation period 2004','FontSize',17,'FontWeight','bold')
ylabel(ax(1),{'atmospheric longwave','radiation [W m^{-2}]'},'FontSize',17,'FontWeight','bold')
ylabel(ax(2),{'above-canopy','air temperature [°C]'},'FontSize',17,'FontWeight','bold')
ylabel(ax(3),{'incoming shortwave','radiation [W m^{-2}]'},'FontSize',17,'FontWeight','bold')
ylabel(ax(4),{'effective emissivity','of the sky'},'FontSize',17,'FontWeight','bold')
xlim(ax(1),[timeticks(1) timeticks(end)]); ylim(ax(1),[-170 470]);
xlim(ax(2),[timeticks(1) timeticks(end)]); ylim(ax(2),[-26 44]);
xlim(ax(3),[timeticks(1) timeticks(end)]); ylim(ax(3),[0 2300]);
xlim(ax(4),[timeticks(1) timeticks(end)]); ylim(ax(4),[-3 1.2]);
set(ax(1),'YTick',150:50:400)
set(ax(2),'YTick',-5:5:20)
set(ax(3),'YTick',0:100:900)
set(ax(4),'YTick',0.6:0.2:1.2)
ylab1 = get(ax(1),'YLabel'); set(ylab1,'Position',get(ylab1,'Position') - [6 -125 0])
ylab2 = get(ax(2),'YLabel'); set(ylab2,'Position',get(ylab2,'Position') + [4 -1.5 0])
ylab3 = get(ax(3),'YLabel'); set(ylab3,'Position',get(ylab3,'Position') - [6 700 0])
ylab4 = get(ax(4),'YLabel'); set(ylab4,'Position',get(ylab4,'Position') + [4 1.78 0])
set(ax(1),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(2),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(3),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(4),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
ax(5) = axes('yaxislocation','left','Color','none','XColor','k','YColor','k');
ax(6) = axes('yaxislocation','right','Color','none','XColor','k','YColor','k');
h5 = line(x,y4,'Parent',ax(5),'Color','k');
h6 = line(x,y5,'Parent',ax(6),'Color','k');
set(ax(5),'Box','off'); set(ax(6),'Box','off');
set(ax(5),'XTick',timeticks,'XTickLabel',datelabel);%set(ax(5),'XTick',[4418 5114 5858 6578]); datetick(ax(5),'x','dd mmm');
set(ax(6),'XTick',timeticks,'XTickLabel',datelabel);%set(ax(6),'XTick',[4418 5114 5858 6578]); datetick(ax(6),'x','dd mmm');
set(ax(5),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(6),'FontSize',17,'FontWeight','bold','LineWidth',2)
xlim(ax(5),[timeticks(1) timeticks(end)]); ylim(ax(5),[0 1]);
xlim(ax(6),[timeticks(1) timeticks(end)]); ylim(ax(6),[0 1]);
set(ax(5),'YTick',[]); set(ax(6),'YTick',[]);
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 20 10];
fig.PaperSize = [20 11];
print(fig,'-dpdf','-r600','Figure3_MetForcEmsky.pdf')

% sub-canopy LW comparison
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
plot([200 450],[200 450],'k')
plot(LW_in_bc_val,LW_in_bc_CLM,'Color',EcstaticEmerald,'Marker','.','LineStyle','none')
plot(LW_in_bc_val,LW_in_bc_SP,'Color',CandidCoral,'Marker','.','LineStyle','none')
xlim([200 450])
ylim([200 450])
xlabel('observed sub-canopy LW radiation [W m^{-2}]','FontSize',17,'FontWeight','bold')
ylabel('simulated sub-canopy LW radiation [W m^{-2}]','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
box on
pbaspect([1 1 1])
fig.PaperUnits = 'inches';
fig.PaperSize = [6.6 6.6];
fig.PaperPosition = [0 0 7 7];
print(fig,'-dpdf','-r600','Figure4_LWbcScatter.pdf')

% effective emissivity dependence
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
scatter(ema_eff,LW_in_bc_CLM_RelErr,17,SW_in_ac_val,'filled','MarkerEdgeColor','none')
colormap(hot(100))
c=colorbar('SouthOutside');
caxis([0 1000])
ylim([-20 20])
xlabel('effective emissivity of the sky','FontSize',17,'FontWeight','bold')
ylabel('sub-canopy LW radiation error [%]','FontSize',17,'FontWeight','bold')
ylabel(c,'incoming SW radiation [W m^{-2}]','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [-0.1 -0.1 8 6];
fig.PaperSize = [7.7 5.9];
set(gcf, 'Renderer', 'opengl')
print(fig,'-dpdf','-r600','Figure5_LWbcScatterEffemaSWac.pdf')

% LW enhancement subplot
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
% a
sp1 = subplot(2,2,1);
hold on
text(0.85,1.65,'a','FontSize',15,'FontWeight','bold')
plot([0.8 1.7],[0.8 1.7],'k')
plot(LWenh_val,LWenh_CLM,'Color',EcstaticEmerald,'Marker','.','Markersize',15,'LineStyle','none')
plot(LWenh_val,LWenh_SP,'Color',CandidCoral,'Marker','.','Markersize',15,'LineStyle','none')
xlim([0.8 1.7])
ylim([0.8 1.7])
% xlabel('observed LW enhancement','FontSize',15,'FontWeight','bold')
ylabel('simulated LW enhancement','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.1:1.7)
set(gca,'YTick',0.8:0.1:1.7)
box on
pbaspect([1 1 1])
% b
sp2 = subplot(2,2,2);
hold on
text(0.85,1.65,'b','FontSize',15,'FontWeight','bold')
plot([0.8 1.7],[0.8 1.7],'k')
plot(LWenh_val,LWenh_CLM_HM,'Color',EcstaticEmerald,'Marker','.','Markersize',15,'LineStyle','none')
%plot(LWenh_val,LWenh_SP,'Color',CandidCoral,'Marker','.','Markersize',15,'LineStyle','none')
xlim([0.8 1.7])
ylim([0.8 1.7])
% xlabel('observed LW enhancement','FontSize',15,'FontWeight','bold')
% ylabel('simulated LW enhancement','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.1:1.7)
set(gca,'YTick',0.8:0.1:1.7)
box on
pbaspect([1 1 1])
% c
sp3 = subplot(2,2,3);
hold on
text(0.85,1.65,'c','FontSize',15,'FontWeight','bold')
plot([0.8 1.7],[0.8 1.7],'k')
scatter(LWenh_val,LWenh_CLM,15,ema_eff,'filled','MarkerEdgeColor','none')
colormap(sp3,x2b2y('single',0.6,1,1.1,'hotncold'))
cb1=colorbar('SouthOutside');
xlim([0.8 1.7])
ylim([0.8 1.7])
xlabel('observed LW enhancement','FontSize',14,'FontWeight','bold')
ylabel('simulated LW enhancement','FontSize',14,'FontWeight','bold')
ylabel(cb1,'effective emissivity of the sky','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.1:1.7)
set(gca,'YTick',0.8:0.1:1.7)
box on
pbaspect([1 1 1])
% d
sp4 = subplot(2,2,4);
hold on
text(0.85,1.65,'d','FontSize',15,'FontWeight','bold')
plot([0.8 1.7],[0.8 1.7],'k')
scatter(LWenh_val,LWenh_CLM,15,SW_in_ac_val,'filled','MarkerEdgeColor','none')
colormap(sp4,hot(100))
cb2=colorbar('SouthOutside');
caxis([0 1000])
xlim([0.8 1.7])
ylim([0.8 1.7])
xlabel('observed LW enhancement','FontSize',14,'FontWeight','bold')
% ylabel('simulated LW enhancement','FontSize',15,'FontWeight','bold')
ylabel(cb2,'incoming SW radiation [W m^{-2}]','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.1:1.7)
set(gca,'YTick',0.8:0.1:1.7)
box on
pbaspect([1 1 1])
% positioning & print
set(sp1,'Position',[0.15,0.625,0.325,0.325])
set(sp2,'Position',[0.575,0.625,0.325,0.325])
set(sp3,'Position',[0.15,0.235,0.325,0.325])
set(sp4,'Position',[0.575,0.235,0.325,0.325])
set(cb1,'Position',[0.15,0.1,0.325,0.05])
set(cb2,'Position',[0.575,0.1,0.325,0.05])
fig.PaperUnits = 'inches';
fig.PaperPosition = [-0.55 -0.1 9 9];
fig.PaperSize = [8 8.8];
set(gcf, 'Renderer', 'opengl')
print(fig,'-dpdf','-r600','Figure6_LWenhScatterSubplotEmptyafter.pdf')

% PDF
spectrum_LWenh = 0.8:0.025:1.6;
spectrum_LWenh_xaxis = 0.8125:0.025:1.5875;
hist_val = histogram(LWenh_val(4225:end),spectrum_LWenh,'Normalization','probability');
    hist_val = hist_val.Values;
hist_CLM = histogram(LWenh_CLM(4225:end),spectrum_LWenh,'Normalization','probability');
    hist_CLM = hist_CLM.Values;
hist_CLM_HM = histogram(LWenh_CLM_HM(4225:end),spectrum_LWenh,'Normalization','probability');
    hist_CLM_HM = hist_CLM_HM.Values;
hist_SP = histogram(LWenh_SP(4225:end),spectrum_LWenh,'Normalization','probability');
    hist_SP = hist_SP.Values;

fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
plot(spectrum_LWenh_xaxis,hist_val,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_CLM,'Color',EcstaticEmerald,'LineWidth',2,'LineStyle',':')
plot(spectrum_LWenh_xaxis,hist_CLM_HM,'Color',EcstaticEmerald,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_SP,'Color',CandidCoral,'LineWidth',2)
hold off
xlim([0.8 1.6])
ylim([0 0.17])
xlabel('LW enhancement','FontSize',17,'FontWeight','bold')
ylabel('Probablity','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.1:1.6)
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 7];
fig.PaperSize = [7.7 6.6];
print(fig,'-dpdf','-r600','Figure7_LWenhPDF.pdf')