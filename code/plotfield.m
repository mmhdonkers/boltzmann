axes('YTick',[0 0.25 0.5 0.75 1],'XTick',[0 0.1 0.2 0.3 0.4 0.5],'Interpreter','latex','FontSize',16)
hold on
axis([0 0.5 0 1])


xlabel('velocity $[m/s]$','Interpreter','latex','FontSize',16)
ylabel('relative tube height','Interpreter','latex','FontSize',16)

a1 = plot(f1,x,'r-');
a2 = plot(f2,x,'b-');
a3 = plot(f3,x,'g-');
a4 = plot(f4,x,'k-');
b1 = plot(vel_0_6_1,x,'ro');
b2 = plot(vel_0_7_1,x,'bo');
b3 = plot(vel_0_8_1,x,'go');
b4 = plot(vel_0_9_1,x,'ko');

set(get(get(b1, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
set(get(get(b2, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
set(get(get(b3, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
set(get(get(b4, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off')
lngd = legend('$\tau=0.6$','$\tau=0.7$','$\tau=0.8$','$\tau=0.9$');
set(lngd,'interpreter','latex','FontSize',16);

hold off