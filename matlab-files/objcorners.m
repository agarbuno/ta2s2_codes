function Hphi = objcorners(phi, a)



Hphi = (4*a - (abs(phi(:,1) - a/2)) - (abs(phi(:,2) - a/2)));
index = 1-prod((phi <= 10 & phi >= 0),2);
if sum(index) > 0
    Hphi(index) = Inf;
end