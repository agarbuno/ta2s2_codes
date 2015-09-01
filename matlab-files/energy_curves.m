function [H, Xz, Yz ] = energy_curves(psi, Nmin, Nmax, nlevels, meanopt, X, y, Temp, H)
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

if isempty(H)
    H = zeros(length(params),1);
    for i=1:length(H)
        H(i) = likelihood(X, y, meanfunc, covfunc, params(i,1:2), ...
            psi, meanopt);
        if imag(H(i)) ~= 0
            H(i) = Inf;
        end
    end
end

L = exp(-H./Temp);
% L = H./Temp;
Lmat = reshape(L, size(Xz));
hold on
%surf(Xz,Yz,Lmat);
[C,h] = contour(Xz,Yz,Lmat,nlevels);
% [C,h] = contour(Xz,Yz,Lmat,linspace(min(H), max(H(H~=Inf)), nlevels));
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('$\theta_1$', 'interpreter', 'latex')
ylabel('$\theta_2$', 'interpreter', 'latex')
title1 = sprintf('Number of training runs = %4i', n);
title2 = sprintf('. Temperature = %4.2f', Temp);
title(strcat(title1, title2),  'interpreter', 'latex')
line([10^Nmin 10^Nmax], [10^Nmin 10^Nmax], 'LineStyle', '--');
hold off
