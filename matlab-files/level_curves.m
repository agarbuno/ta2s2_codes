function level_curves(psi, Nmin, Nmax, nlevels, meanopt, X, y, lb)

global n;

% Gaussian process settings
meanfunc = 'h';
covfunc = 'covK';

% Prior distribution for gaussian parameters
prior = 'randn';

% Grid size
N = 70;

% Generate test points in log scale
theta1 = logspace(Nmin,Nmax,N)';
theta2 = logspace(Nmin,Nmax,N)';

[Xz, Yz] = meshgrid(theta1, theta2);
params = [ reshape(Xz,N^2,1) reshape(Yz,N^2,1) ];

L = zeros(length(params), 1);

for i=1:length(L)
    L(i) = likelihood(X, y, meanfunc, covfunc, params(i,1:2), ...
        psi, meanopt, lb);
    if imag(L(i) )
        flag = 0;
        L(i) = Inf;
    end
    if mod(i,10000) == 0
        fprintf(1,'Iteration %12i \n', i);
    end
end

fprintf(1,'\nMinimum empiric likelihood level: %4.2f\n\n', min(L));

Lmat = reshape(L, size(Xz));
% Lmat = exp(-Lmat);

hold on
% surf(Xz,Yz,Lmat);
contour3(Xz,Yz,Lmat,nlevels);
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('$\phi_1$', 'interpreter', 'latex')
ylabel('$\phi_2$', 'interpreter', 'latex')
title(sprintf('Number of training runs = %4i', n),  'interpreter', 'latex')
% line([10^Nmin 10^Nmax], [10^Nmin 10^Nmax], 'LineStyle', '--');
hold off
% pause
% Printing options
% set(gcf, 'PaperPosition', [0.6350 6.3500 20.3200 15.2400]);
% set(gcf, 'PaperSize', [21.0000 29.7000]); 
% title = sprintf('~/Dropbox/Phd/Research/Reports/Integrated_likelihood/level_curves_sample%03d.eps' ,n);
% saveas( gcf, title, 'eps2c')

% hold on
% surf(Xz,Yz,Lmat);
% hold off
% 
% figure(2); clf;
% surf(Xz,Yz,log(log(Lmat)));
% set(gca,'yscale','log');
% set(gca,'xscale','log');
% xlabel('$\theta_1$', 'interpreter', 'latex')
% ylabel('$\theta_2$', 'interpreter', 'latex')

% contour(log(Xz),log(Yz),Lmat,70);

% covfunc = 'covKiso';
% gamma = .5;
% alpha = 0.05;
% N = 1000;
% c = 1;
% meanopt = 3;
% 
% tic
%     [beta, theta, sigma, Hnew, k, w] = aims_opt(X, y, gamma, alpha, N, c, ...le
%         meanfunc, covfunc, meanopt, prior);
% toc
% 
% fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));
% fprintf(1, 'Minimum : [ %4.8f ] \n', min(L));
% 
% scatter(exp(theta), exp(sigma))
% hold off