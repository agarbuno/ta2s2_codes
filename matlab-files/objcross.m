function Hphi = objcross(phi, a)

Hphi = 1 + ((phi(:,1) - a/2).^2) .* ((phi(:,2) - a/2).^2 );