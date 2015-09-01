function Hphi = objcrosssamp(phi, a)

N = 1e3;
Hphi = zeros(size(phi,1),1);

for i = 1:size(phi,1)
    theta = bsxfun(@plus, randn(N,2), (phi(i,:) - 0.5 * a));
    Hphi(i) = 1 + mean(prod(phi(i,:) - a/2, 2) .* prod(theta,2));
end

% Hphi = 1 + ((phi(:,1) - a/2).^2) .* ((phi(:,2) - a/2).^2 )