function output = Laguerre(X,k) % k - order of laguerre, X - input vector
    if k==0
        output = exp(-X/2);
    elseif k == 1
        output = exp(-X/2).*(1-X);
    elseif k == 2
        output = exp(-X/2).*(1-2*X+X.^2/2);
    end
end