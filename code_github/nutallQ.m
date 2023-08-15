%
%   Q_{n+2k+1,n}(a,b) 
%
function out = nutallQ(m, n, a, b)
    %
    k = (m-n-1)/2;
    if (isinteger(k))
        %
        if k >= 2
            out = a*nutallQ(n+2*k, n+1, a, b) + 2*(n+k)*nutallQ(n+2*k-1, n, a, b)...
                + b^(n+2*k) * exp( -a^2/2-b^2/2 ) * besseli(n, a*b);
        else % k = 1
            out = a^(n+2) * marcumq(a, b, n+2) + 2*(n+1)*a^n*marcumq(a, b, n+1)...
                + b^(n+2) * exp( -a^2/2-b^2/2 ) * besseli(n, a*b);
        end
    else
        Infinity = 50; % 
        fx = @(x) x.^m .* exp(-x.^2/2-a^2/2).*besseli(n, a*x);
        out = integral(fx, b, Infinity);
    end
end