function h = h(x,opt)

[n,m] = size(x);

switch opt
    case 1
        h = ones(n,1);
    case 2
        h = [ones(n,1), x];
    case 3
        h = zeros(n,1); 
    case 4
        h = ones(n,1);
        for i = 1:2
           h = [ h, x.^i];
        end
    otherwise
        h = [];
        fprintf(1,'Not a valid option for h(x)\n')
end


