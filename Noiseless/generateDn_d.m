% construction from Hayashi and Ouyang
function D = generateDn_d(n, d)
    if d == 1
        theta = 2*pi*(0:n-1)'/n;
        D = [cos(theta), sin(theta)];   
        return;
    end
    D_prev = generateDn_d(n, d-1);  
    D = zeros(n * size(D_prev,1), d+1);
    idx = 1;
    for j = 0:n-1
        alpha = 2*pi*j/n;
        c = cos(alpha);
        s = sin(alpha);
        block = [c * D_prev, s * ones(size(D_prev,1), 1)];
        D(idx : idx + size(D_prev,1) - 1, :) = block;
        idx = idx + size(D_prev,1);
    end
end