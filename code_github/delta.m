function out = delta(n, a, b, s)

    c = b./(1+a);
    out = factorial(n-1)*(-c./(1+c*s)).^n.*(n*a./(1+c*s)+1);
end