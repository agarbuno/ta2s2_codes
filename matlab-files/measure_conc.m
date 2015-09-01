close all
Nsamp = 1000;
figure(1); clf; hold on
for i = 0:4
    d = 10^i; x = randn( Nsamp, d); 
    h = histogram(sqrt(sum(x.^2,2))/sqrt(d));     
    h.Normalization = 'pdf';
    h.Parent.YScale = 'log';
    h.FaceAlpha = (1+i)/5;
    h.FaceColor = 'blue';
    pause
end
hold off

close all
Nsamp = 5000;
figure(1); clf; hold on
for i = 0:4
    d = 10^i; x = rand( Nsamp, d); 
    index = randsample(1:Nsamp,1);
    point = x(index,:); x = x(1:Nsamp ~= index, :);
    h = histogram(sqrt(sq_dist(point', x'))/sqrt(d));
    h.Normalization = 'pdf';
    h.Parent.YScale = 'log';
    h.FaceAlpha = (1+i)/5;
    h.FaceColor = 'red';
    pause
end
hold off

close all
Nsamp = 10000;
figure(1); clf; hold on
for i = 0:4
    
    d = 10^i; x = -1/2 + rand( Nsamp, d); 
    figure(1); 
    h =  histogram( (sqrt(sum(abs(x),2))/( sqrt(d)/2 )) , 20);  
    h.Normalization = 'pdf';
    h.Parent.YScale = 'log';
    h.FaceAlpha = (1+i)/5;
    h.FaceColor = [255/255;127/255;36/255];
     
    pause
end
hold off
