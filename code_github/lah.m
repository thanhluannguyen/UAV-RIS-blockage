function out = lah(n,k)
    if (n == 0) && (k == 0)
        c = 1;
    else
        c = k/n;
    end
    out = (factorial(n) ./ factorial(k)).^2 .* c ./factorial(n-k);
end