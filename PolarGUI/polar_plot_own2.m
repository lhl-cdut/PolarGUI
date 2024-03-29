function [] = polar_plot_own2(Pazi4)
% Pazi4:Incidence 0-180

%%
[r,theta] = histcounts(Pazi4);
theta = theta + (theta(2)-theta(1))/2;
r_max = max(r);
%%
% subplot(2,2,1);

for i=1:1:length(r)
    px = [0 r(i)* sind(theta(i))];
    py = [0 r(i)* cosd(theta(i))];
    
    plot(px,py,'-','linewidth',2,'Color','[0.00,0.35,0.67]');
    clear px py
    hold on;
end

cNum = r_max/3;
if cNum > 2
    for i=0:cNum:r_max
        text(i* sind(10),i* cosd(10),num2str(ceil(i)),'Color','k');
        plot(i* sind(0:360),i* cosd(0:360),'--','Color',[0.80,0.80,0.80])
        hold on
    end
else
    for i=0:2:r_max
        text(i* sind(10),i* cosd(10),num2str(ceil(i)),'Color','k');
        plot(i* sind(0:360),i* cosd(0:360),'--','Color',[0.80,0.80,0.80]);
        hold on
    end
end

plot(r_max* sind(0:360),r_max* cosd(0:360),'Color',[0.80,0.80,0.80]);
hold on;

for i=0:30:330
    plot([0 r_max* sind(i)],[0 r_max* cosd(i)],'--','Color',[0.80,0.80,0.80]);
    hold on;
    
    if i==0 || i==90
       continue; 
    end
    
    if i<180
        text((r_max+cNum/2)* sind(i),(r_max+cNum/2)* cosd(i),[num2str(i),'^��'],'Color','k');
    elseif i==180
        text(r_max*sind(190),r_max*cosd(190)-cNum/2,[num2str(i),'^��'],'Color','k');
    elseif i>180 && i<270
        text(r_max* sind(i)-cNum*1.5,r_max* cosd(i)-cNum/2,[num2str(i),'^��'],'Color','k');
    elseif i==270
        text(r_max* sind(i)-cNum*1.5,r_max* cosd(i),[num2str(i),'^��'],'Color','k');
    else
        text(r_max* sind(i)-cNum*1.5,r_max* cosd(i)+cNum/2,[num2str(i),'^��'],'Color','k');
    end
end
text(r_max*sind(345),(r_max+cNum/1.5),'North','Color','k');
text(r_max+cNum/5,0,'East','Color','k');

% grid on;
% title(titlename,'color','w');
set(gca,'XLim',[-r_max-cNum*0.5 r_max+cNum*0.5],'YLim',[-r_max-cNum*0.5 r_max+cNum*0.5]);
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
% set(gca,'DataAspectRatio',[1 1 1]);
axis equal;

end