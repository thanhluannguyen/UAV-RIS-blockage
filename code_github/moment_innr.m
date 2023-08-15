function out = moment_innr(M, K_aerial, inr_aerial, K_ground, inr_ground, x)
    %
    s = 0;
    eps = 1;
    %
    syms p
    laplace_aerial = prod( (1+K_aerial)./(1+K_aerial+inr_aerial*p)...
        .* exp(-K_aerial.*inr_aerial*p./(1+K_aerial+inr_aerial*p)) );
    laplace_ground = prod( (1+K_ground)./(1+K_ground+inr_ground*p)...
        .* exp(-K_ground.*inr_ground*p./(1+K_ground+inr_ground*p)) );

    laplace_cci_awgn_p = exp(-p) * laplace_aerial * laplace_ground;
    laplace_cci_awgn_s = double(limit(laplace_cci_awgn_p, p, s));
    
    SUM_i = [];
    for n = 1:M
        temp = sum(delta(n, K_aerial, inr_aerial, p))...
             + sum(delta(n, K_ground, inr_ground, p));
        if n == 1
            SUM_i = [SUM_i, (-eps*x)^n * (temp-1)];
        else
            SUM_i = [SUM_i, (-eps*x)^n * (temp)];
        end
    end
    SUM_i = double(limit(SUM_i, p, s));
    % loading the record (see genRecord.m)
    load('record.mat');
    %
    sum_n = 0; upper = size(record, 1);
    %
    n = min(M, upper); % this algorithm cannot run for all values of M
    %
    diff_lap_cci_awgn = sum( record{n, 1}.*prod((SUM_i(record{n, 2})).^record{n, 3}, 1), 'all' );
    %
    out = laplace_cci_awgn_s * diff_lap_cci_awgn / eps^n;
end