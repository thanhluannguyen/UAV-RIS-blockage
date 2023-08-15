function [param, PDF, CDF] = fitdist_SRS_asymp(mu_1, mu_2, m, Q)
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
    f_asymp = @(x) (alp/(xi+1))^(Q*m-Q)/factorial(Q-1) * alp^(Q*m) * x.^(Q-1) .* exp(-x/alp);
    %
    %
    F_asymp = @(x) 0;
    for k = 0:(Q*m-Q)
        chi_Nk = nchoosek(Q*m-Q, k)/(xi+1)^k * (xi/(xi+1))^(Q*m-Q-k);

        F_asymp = @(x) F_asymp(x) + chi_Nk .* (x/alp).^(Q*m-k)/factorial(Q*m-k);
    end
    %
%     F_asymp = @(x) xi^(-Q*m) * (xi/(xi+1))^(Q*m-Q) * (x*xi/alp).^Q / factorial(Q);
    %
    param = [kappa omega xi Omega alp Q];
    PDF = @(x) f_asymp(x);
    CDF = @(x) F_asymp(x);
    %
    % figure;
    % fplot(CDF, [0, max(Y)]); hold off; 
    % set(gca, 'XScale', 'Log');
    % set(gca, 'YScale', 'Log');
    %
end