load('solvedData');
load('fitPeaks');
numScans = 23;

yaxlims = [-0.25,0.62;-0.25 0.25;-.03 .5;-2.6 1.1];
totyrange = sum(yaxlims(:,2)-yaxlims(:,1));
figgap = 0.03;
availplotrange = 0.95 - 0.06 - figgap * 3;
figlocsy = 0.06 + [0; cumsum((yaxlims(:,2)-yaxlims(:,1))/totyrange * availplotrange)] ...
    + [0; figgap; 2*figgap; 3*figgap; 4*figgap];

posax = (posax + posax(end-1) - 2* posax(end))*1000;

colorvec = sunset;
colors = [
    1.0000    0.4980    0.0157
    0.7451    0.2824    0.4196
    0.9294    0.3804    0.2824
    0.4431    0.1882    0.4549
    0.1412    0.0941    0.4863
    0.3843    0.4078    0.6980];

solPeaksPlot = reshape(solPeaks,[16,numScans]);
peaksPlot = reshape(centers,[16,numScans]) * 1000;

sol = sol * 10^5;
ranges = ranges * 10^5;

epsList = {'\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}','\epsilon_{xy}','\epsilon_{yz}','\epsilon_{zx}'};

figure(50);clf;
set(gcf,'units','centimeters','position',[0,10,8.2,14.15]);
subplot('Position',[0.15,figlocsy(4),0.8,figlocsy(5)-figlocsy(4)]);
mults = 1:3;
for i = mults
    errorbar(posax(:),sol(i,:),sol(i,:) - squeeze(ranges(i,1,:))',squeeze(ranges(i,2,:))' - sol(i,:),'.','Color',colors(i,:),'CapSize',2.5,'LineWidth',0.8,'MarkerSize',8);
    hold on;
end

set(gca,'fontsize',12);
box off;
set(gca,'XTick',[0 500 1000],'Linewidth',1.5,'XTickLabel',[],'TickDir','out');
ylim([-2.6 1.1]);
ylbl = ylabel('Strain \times10^{-5} (unitless)');
ylbl.Position = [-180,-2.2,0];

for single = 4:6
    if single == 4
        subplot('Position',[0.15,figlocsy(3),0.8,figlocsy(4)-figlocsy(3)-figgap]);
    elseif single == 5
        subplot('Position',[0.15,figlocsy(2),0.8,figlocsy(3)-figlocsy(2)-figgap]);
    else
        subplot('Position',[0.15,figlocsy(1),0.8,figlocsy(2)-figlocsy(1)-figgap]);
    end
    if single == 4
        errorbar(posax(:),...
            abs(sol(single,:)),...
            abs(sol(single,:)),...
            max(squeeze(ranges(single,2,:))',-squeeze(ranges(single,1,:))') - abs(sol(single,:)),...
            '.','Color',colors(single,:),'CapSize',2.5,'LineWidth',0.8,'MarkerSize',8);
    else
        errorbar(posax(:),...
            sol(single,:),...
            sol(single,:) - squeeze(ranges(single,1,:))',...
            squeeze(ranges(single,2,:))' - sol(single,:),...
            '.','Color',colors(single,:),'CapSize',2.5,'LineWidth',0.8,'MarkerSize',8);
    end
    box off
    if single == 4
        ylim([-.03 .5])
        set(gca,'XTick',[0 500 1000],'Linewidth',1.5,'XTickLabel',[],'TickDir','out');
        set(gca,'YTick',[0 .25 .5]);
    elseif single == 5
        ylim([-0.25 0.25]);
        set(gca,'XTick',[0 500 1000],'Linewidth',1.5,'XTickLabel',[],'TickDir','out');
        set(gca,'YTick',[-.25 0 .25]);
    else
        ylim([-0.25,0.62]);
        set(gca,'XTick',[0 500 1000],'Linewidth',1.5,'TickDir','out');
        set(gca,'YTick',[-.25 0 .25 .5]);
    end
end
xlabel('Position (μm)');

set(findall(gcf,'-property','FontSize'),'FontSize',8)
set(gcf, 'Color', 'w');