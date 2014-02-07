axes('YTick',[0 0.25 0.5 0.75 1],'XTick',[0 0.05 0.1 0.15 0.2])
hold on
axis([0 0.2 0 1])


xlabel('velocity $[m/s]$','Interpreter','latex','FontSize',14)
ylabel('relative tube height','Interpreter','latex','FontSize',14)

ra1 = plot(rf1,x,'r-');
ra2 = plot(rf2,x,'b-');
ra3 = plot(rf3,x,'g-');
ra4 = plot(rf4,x,'k-');
ra5 = plot(rf5,x,'y-');
rb1 = plot(vel_0_9_0_6,x,'ro');
rb2 = plot(vel_0_9_0_8,x,'bo');
rb3 = plot(vel_0_9_1,x,'go');
rb4 = plot(vel_0_9_1_2,x,'ko');
rb5 = plot(vel_0_9_1_4,x,'yo');

set(get(get(rb1, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
set(get(get(rb2, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
set(get(get(rb3, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
set(get(get(rb4, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
set(get(get(rb5, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
lngd = legend('$\rho=0.6$','$\rho=0.8$','$\rho=1$','$\rho=1.2$','$\rho=1.4$');
set(lngd,'interpreter','latex','FontSize',14,'Location','West');

hold off