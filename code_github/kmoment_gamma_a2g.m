function out = kmoment_gamma_a2g(k, mu_cascaded, mu_direct)

    temp = 0;
    for t = 0:(2*k)
        if (t == 2*k), A = 1;
        else
            A = mu_cascaded(2*k-t);
        end
        if (t == 0), B = 1;
        else
            B = mu_direct(t);
        end
        temp = temp + nchoosek(2*k, t) * A * B;
    end

    out = temp;
end