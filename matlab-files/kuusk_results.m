% load('results_kuusk_paper')
figure(1); clf; hold on; ksdensity(Hnew, 'width', 0.25); set(gca, 'FontSize', 15); 
xlabel('$\mathcal{H}(\phi)$', 'interpreter', 'latex')
ylabel('Estimated pdf', 'interpreter', 'latex')
title('Estimated density of $\mathcal{H}(\phi)$',  'interpreter', 'latex')
hold off

i = 1;
figure(i+1); ksdensity(exp(theta(:,i))); set(gca, 'FontSize', 15); 
xlabel('$\phi_1$', 'interpreter', 'latex')
ylabel('Estimated pdf', 'interpreter', 'latex')
title('Estimated density of $\phi_1$',  'interpreter', 'latex')

i = 2;
figure(i+1); ksdensity(exp(theta(:,i))); set(gca, 'FontSize', 15); 
xlabel('$\phi_2$', 'interpreter', 'latex')
ylabel('Estimated pdf', 'interpreter', 'latex')
title('Estimated density of $\phi_2$',  'interpreter', 'latex')

i = 3;
figure(i+1); ksdensity(theta(:,i), 'width', 0.1); set(gca, 'FontSize', 15); 
xlabel('$\log(\phi_3)$', 'interpreter', 'latex')
ylabel('Estimated pdf', 'interpreter', 'latex')
title('Estimated density of $\phi_3$',  'interpreter', 'latex')

i = 4;
figure(i+1); ksdensity(theta(:,i), 'width', .2); set(gca, 'FontSize', 15); 
xlabel('$\log(\phi_4)$', 'interpreter', 'latex')
ylabel('Estimated pdf', 'interpreter', 'latex')
title('Estimated density of $\phi_4$',  'interpreter', 'latex')

i = 5;
figure(i+1); ksdensity(exp(theta(:,i)), 'width', .0007); set(gca, 'FontSize', 15); 
xlabel('$\phi_5$', 'interpreter', 'latex')
ylabel('Estimated pdf', 'interpreter', 'latex')
title('Estimated density of $\phi_5$',  'interpreter', 'latex')
