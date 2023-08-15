function out = cdf_sinr(Q, m, kappa, omega, K_aerial, inr_aerial, K_ground, inr_ground, x)
    %
        xi = kappa*omega/m;
        alp = (xi+1)/(kappa+1);
    s = x / alp;

    syms p
    laplace_aerial = prod( (1+K_aerial)./(1+K_aerial+inr_aerial*p)...
        .* exp(-K_aerial.*inr_aerial*p./(1+K_aerial+inr_aerial*p)) );
    laplace_ground = prod( (1+K_ground)./(1+K_ground+inr_ground*p)...
        .* exp(-K_ground.*inr_ground*p./(1+K_ground+inr_ground*p)) );

    laplace_cci_awgn_p = exp(-p) * laplace_aerial * laplace_ground;
    laplace_cci_awgn_s = double(limit(laplace_cci_awgn_p, p, s));
    
    S_i = [];
    for n = 1:Q*m-1
        temp = sum(delta(n, K_aerial, inr_aerial, p))...
             + sum(delta(n, K_ground, inr_ground, p));
        if n == 1
            S_i = [S_i (-p)^n * (temp-1)];
        else
            S_i = [S_i (-p)^n * (temp)];
        end
    end
    S_i = double(limit(S_i, p, s));
    % loading the record (see genRecord.m)
    load('record.mat');
    %
    data_diff_cci_awgn = NaN * ones(Q*m-1, 1);
    %
    sum_k =  0;
    for k = 0:Q*m-Q
        chi = nchoosek(Q*m-Q, k)/(xi+1)^k * (xi/(xi+1))^(Q*m-Q-k);

        sum_n = 0; upper = size(record, 1);
        for n = 0:min(Q*m-k-1, upper)
            if n > 0
                if isnan(data_diff_cci_awgn(n))
                    diff_lap_cci_awgn = sum( record{n, 1}.*prod((S_i(record{n, 2})).^record{n, 3}, 1), 'all' );
                    data_diff_cci_awgn(n) = diff_lap_cci_awgn;
                else
                    diff_lap_cci_awgn = data_diff_cci_awgn(n); 
                end
            else
                diff_lap_cci_awgn = 1;
            end

            % SUM_i(1)^2 + SUM_i(2)
            % SUM_i(1)^3 + 3*SUM_i(2)*SUM_i(1) + SUM_i(3)
            % SUM_i(1)^4 + 6*SUM_i(2)*SUM_i(1)^2 + 4*SUM_i(1)*SUM_i(3) + 3*SUM_i(2)^2 + SUM_i(4)
            % SUM_i(1)^5 + 10*SUM_i(1)^3*SUM_i(2) + 10*SUM_i(1)^2*SUM_i(3) ... 
            % + 15*SUM_i(1)*SUM_i(2)^2 + 5*SUM_i(1)*SUM_i(4) + 10*SUM_i(2)*SUM_i(3) + SUM_i(5)

            sum_n = sum_n + laplace_cci_awgn_s * diff_lap_cci_awgn;
        end
        sum_k =  sum_k + chi * sum_n;
    end

    out =  1 - sum_k;
end