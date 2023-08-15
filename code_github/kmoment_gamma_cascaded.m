%
% Calculate the kth moment of sqrt(gamma_cascaded)
%
function out = kmoment_gamma_cascaded(N, k, K1, K2)
    %
    mu = @(k) (1+K1).^(-k/2).*gamma(1+k/2).*laguerreL(k/2, -K1)...
           .* (1+K2).^(-k/2).*gamma(1+k/2).*laguerreL(k/2, -K2);
    %
    M = generate_partitions(k);
    %
    %  M contains all possible integer partitions of (k), 
    %   where len(M{i}) = r
    %
    sum_m = 0;
    for iM = 1:length(M)
        [eigs, ms] = unique_eigenvalues(diag(M{iM}));
        
        r = length(M{iM});

        sum_m = sum_m + prod(mu(eigs).^ms./factorial(ms)) / prod(factorial(M{iM}))...
            * nchoosek(N, r)*factorial(r) * factorial(k);
    end

    out = sum_m;
end