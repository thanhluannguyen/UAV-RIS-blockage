function cdf_asy = cdf_asy_sinr(mu_1, mu_2, m, Q, K_aerial, inr_aerial, K_ground, inr_ground, innr, x)
    %
    warning off
    %
    max_N = floor(mu_1^2/(mu_2-mu_1^2));
    %
    if Q > max_N
        fprintf('Q should not exceed: %d \n', max_N);
    end 
    Q = min(Q, max_N);
    %
    Omega = sqrt(m/(m-1)) * sqrt(mu_1^2*(Q+1)-Q*mu_2)/Q;
    sigma2= mu_1/Q - Omega;
    kappa = 1/sigma2-1;
    omega = Omega * (kappa+1)/kappa;
    xi = kappa*omega/m;
    alp = (xi+1)/(kappa+1);
    %
    F_asymp =  0;
    for k = 0:(Q*m-Q)
        chi_Nk = nchoosek(Q*m-Q, k)/(xi+1)^k * (xi/(xi+1))^(Q*m-k-Q);
        
%         F_asymp = F_asymp + chi_Nk...
%             .* mean( (x*innr/alp).^(Q*m-k) ) / factorial(Q*m-k);
        %
        F_asymp = F_asymp + chi_Nk...
            .* moment_innr(Q*m-k, K_aerial, inr_aerial, K_ground, inr_ground, x/alp);
    end
%     F_asymp = 1/(xi+1)^(Q*m-Q) * moment_innr(Q, K_aerial, inr_aerial, K_ground, inr_ground, x/(alp*snr));
    %
    cdf_asy = F_asymp;
    %
end