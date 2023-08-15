clear all;
trials = 5e6;
%%
L = 1;      % number of aerial CCI
K = 1;      % number of ground CCI
M = 4;      % number of transmit antenna
N = 36;     % number of reflecting elements
%% Noise
B = 10*10^6; % Hz - Bandwidth
N0= -174;    % dBm/Hz - Noise Spectral Density
AWGN = db2pow(N0)*B; % mW - Noise Power
%% 3D deployment
max_d = 1;
%
xS = 0;
yS = 0;
zS = 0;

% UAV
x_US = 0.5*max_d;
y_US = 0.5*max_d;
z_US = 1*max_d;

% DES
x_DS = 0.5*max_d;
y_DS = 0.5*max_d;
z_DS = 0;

% RIS
x_RS = 1.0*max_d;
y_RS = 1.0*max_d;
z_RS = 0*max_d;

% ACCI
if (L == 2)
    x_IS= [0.2; 0.4]*max_d;
    y_IS= [0.2; 0.8]*max_d;
    z_IS= [0.6; 0.4]*max_d;
else
    x_IS= 0.2*max_d;
    y_IS= 0.2*max_d;
    z_IS= 0.6*max_d;
end

% GCCI
if K == 2
    x_JS= [0.4; 0.9]*max_d;
    y_JS= [0.4; 0.9]*max_d;
    z_JS= [0; 0]*max_d;
else
    x_JS= 0.4*max_d;
    y_JS= 0.4*max_d;
    z_JS= 0*max_d;
end
%% Distance
d_SU = sqrt((xS-x_US).^2 + (yS-y_US).^2 + (zS-z_US).^2);
d_UR = sqrt((x_RS-x_US).^2 + (y_RS-y_US).^2 + (z_RS-z_US).^2);
d_UD = sqrt((x_DS-x_US).^2 + (y_DS-y_US).^2 + (z_DS-z_US).^2);
d_IU = sqrt((x_IS-x_US).^2 + (y_IS-y_US).^2 + (z_IS-z_US).^2);
d_JU = sqrt((x_JS-x_US).^2 + (y_JS-y_US).^2 + (z_JS-z_US).^2);
d_ID = sqrt((x_IS-x_DS).^2 + (y_IS-y_DS).^2 + (z_IS-z_DS).^2);
d_JD = sqrt((x_JS-x_DS).^2 + (y_JS-y_DS).^2 + (z_JS-z_DS).^2);
d_RD = sqrt((x_RS-x_DS).^2 + (y_RS-y_DS).^2 + (z_RS-z_DS).^2);
%% Large-Scale Fading (Pathloss)
Gt = db2pow(0); % 0 - dBi
Gr = db2pow(0); % 0 - dBi
fc = 3; % GHz
PL = @(d) pathloss3GPP_UMi(Gr,Gt,fc,d);
PL_SU = PL(d_SU);
PL_UR = PL(d_UR);
PL_UD = PL(d_UD);
PL_IU = PL(d_IU);
PL_JU = PL(d_JU);
PL_ID = PL(d_ID);
PL_JD = PL(d_JD);
PL_RD = db2pow(- Gr - Gt - 37.3 + 26*log10(fc) - 36.7*log10(d_RD));
%% Small-Scale Fading Constants
K0 = 10^(0/10);
Kpi= 10^(5/10);
K1 = K0;
K2 = 2/pi*log(Kpi/K0);

kappa_ab = @(za,zb,dab) asin(abs(za-zb)./dab); % works for elevation in [-pi/2, pi/2]
kappa_SU = kappa_ab(zS,z_US,d_SU);
kappa_UR = kappa_ab(z_RS, z_US, d_UR);
kappa_UD = kappa_ab(z_DS, z_US, d_UD);
kappa_IU = kappa_ab(z_IS, z_US, d_IU);
kappa_JU = kappa_ab(z_JS, z_US, d_JU);
kappa_ID = kappa_ab(z_IS, z_DS, d_ID);
kappa_JD = kappa_ab(z_JS, z_DS, d_JD);
kappa_RD = kappa_ab(z_RS, z_DS, d_RD);

K_SU = K1*exp(K2*kappa_SU); % Rician factor of gamma_SU
K_UR = K1*exp(K2*kappa_UR); % Rician factor of gamma_UR
K_UD = K1*exp(K2*kappa_UD); % Rician factor of gamma_UD
K_IU = K1*exp(K2*kappa_IU); % Rician factor of gamma_IU
K_JU = K1*exp(K2*kappa_JU); % Rician factor of gamma_JU
K_ID = K1*exp(K2*kappa_ID); % Rician factor of gamma_ID
K_JD = K1*exp(K2*kappa_JD); % Rician factor of gamma_JD
K_RD = K1*exp(K2*kappa_RD); % Rician factor of gamma_RD

lambda_SU = (K_SU + 1)/PL_SU;
lambda_UR = (K_UR + 1)/PL_UR;
lambda_UD = (K_UD + 1)/PL_UD;
lambda_IU = (K_IU + 1)/PL_IU;
lambda_JU = (K_JU + 1)/PL_JU;
lambda_ID = (K_ID + 1)/PL_ID;
lambda_JD = (K_JD + 1)/PL_JD;
lambda_RD = (K_RD + 1)/PL_RD;
%% Normalized Cascaded Gain & Normalized SNR (2st hop)
filename = sprintf('gamma_cascaded_N%d.mat', N);
try
    load(filename);
catch
    gamma_cascaded = 0;
    for n = 1:N
        %
        gamma_UR = ncx2rnd(2, 2*K_UR, [1,trials])/(2*(K_UR+1));
        gamma_RD = ncx2rnd(2, 2*K_RD, [1,trials])/(2*(K_RD+1));
        %
        gamma_cascaded = gamma_cascaded + sqrt(gamma_UR.*gamma_RD);
    end
    gamma_cascaded = gamma_cascaded.^2;
    %
    save(filename, 'gamma_cascaded');
end
%
%% Normalized A2G SNR
ibeta = 1;

generate_marcumqTable(K_UD);
marcumqTable = cell2mat(struct2cell(load('marcumqTable.mat')));
inx = (1:9:100);
marcumqTable = marcumqTable(inx, :);
betas = marcumqTable(:, 1);
% 0
% 0.0100670473331284
% 0.0597571572880559
% 0.191426537753562
% 0.416100628846570
% 0.669431828304073
% 0.860516311805993
% 0.957572629260874
% 0.990891287023915
% 0.998638110850211
% 0.999859415115188
% 0.999990038049460

beta = betas(ibeta);
tau = marcumqTable(ibeta, 2)^2 / (2*lambda_UD);

fprintf('beta = %f, tau = %f dB \n', beta, pow2db(tau));    

filename = sprintf('gamma_UD_%d.mat', ibeta);
try
    load(filename);
catch
    
end

rlfCoeff = db2pow(0); 
bg_cascaded  = (rlfCoeff^2*PL_UR*PL_RD);
gamma_t = PL_UD*gamma_UD/bg_cascaded;

gamma_a2g = sqrt(gamma_cascaded) + sqrt(gamma_t);
gamma_a2g = gamma_a2g.^2;
%
%% Normalized G2A SNR
gamma_g2a = sum(ncx2rnd(2, 2*K_SU, [M, trials])/(2*(K_SU+1)), 1);
%% Aerial Interference Power
PI = db2pow(0*ones(L,1));
bg_IU = PI/AWGN.*PL_IU;
bg_ID = PI/AWGN.*PL_ID;
gamma_IU = 0;
gamma_ID = 0;
for ii = 1:L
    %
    gamma_IU = gamma_IU + bg_IU(ii)*ncx2rnd(2, 2*K_IU(ii), [1, trials])/(2*(K_IU(ii)+1));
    gamma_ID = gamma_ID + bg_ID(ii)*ncx2rnd(2, 2*K_ID(ii), [1, trials])/(2*(K_ID(ii)+1));
end
%% Ground Interference Power
PJ = db2pow(0*ones(K,1));
bg_JU = PJ/AWGN.*PL_JU;
bg_JD = PJ/AWGN.*PL_JD;
gamma_JU = 0;
gamma_JD = 0;
for ii = 1:K
    %
    gamma_JU = gamma_JU + bg_JU(ii)*ncx2rnd(2, 2*K_JU(ii), [1, trials])/(2*(K_JU(ii)+1));
    gamma_JD = gamma_JD + bg_JD(ii)*ncx2rnd(2, 2*K_JD(ii), [1, trials])/(2*(K_JD(ii)+1));
end
%% Normalized SINR G2A
bg_g2a = 1/AWGN*PL_SU;
sinr_g2a = bg_g2a*gamma_g2a./(gamma_IU + gamma_JU + 1);
%% Normalized SINR A2G
sinr_a2g = 1/AWGN*bg_cascaded*gamma_a2g./(gamma_ID + gamma_JD + 1);
sinr_a2g_asy = 1/AWGN*bg_cascaded*mean(gamma_a2g)./( mean(gamma_ID) + mean(gamma_JD) + 1);
%% PDF & CDF of gamma_SU
m_g2a = 4;
Nmax_g2a = M;

mu_gamma_XY = @(k, K) (1+K).^(-k).*gamma(1+k).*laguerreL(k, -K);
mu_g2a_1 = M*mu_gamma_XY(1, K_SU);
mu_g2a_2 = M*mu_gamma_XY(2, K_SU) + (M-1)/M*mu_g2a_1^2;

disp('Matching G2A SNR to SSR:...')
[param_g2a, pdf_gamma_g2a, cdf_gamma_g2a] = fitdist_SRS(mu_g2a_1, mu_g2a_2, m_g2a, Nmax_g2a);
[~, ~, cdf_asy_gamma_g2a] = fitdist_SRS_asymp(mu_g2a_1, mu_g2a_2, m_g2a, Nmax_g2a);
%% PDF & CDF of G2A SINR
kappa_g2a = param_g2a(1);
omega_g2a = param_g2a(2);
xi_g2a = param_g2a(3);
alp_g2a = param_g2a(5);
N_g2a = param_g2a(6);

disp('Calculating G2A SINR...')
cdf_sinr_g2a = @(x) cdf_sinr(N_g2a, m_g2a, kappa_g2a, omega_g2a, K_IU, bg_IU, K_JU, bg_JU, x/bg_g2a);
cdf_asymp_sinr_g2a = @(x) cdf_asy_sinr(mu_g2a_1, mu_g2a_2, m_g2a, N_g2a, K_IU, bg_IU, K_JU, bg_JU, gamma_IU + gamma_JU + 1, x/bg_g2a);
%% PDF & CDF of A2G SNR
f_gamma_UD = @(x) exp(-lambda_UD*x-K_UD).*lambda_UD.*besseli(0, 2*sqrt(K_UD*lambda_UD*x));
Infinity = 2*1e6;
F_gamma_UD = @(x) 1 - marcumq(sqrt(2*K_UD), sqrt(2*lambda_UD*x));

% Transition Matrix
P = [beta 1-beta; beta 1-beta];
% Steady-State Probabilities
pi1 = P(1,2)/(P(1,2)+P(2,1));
pi0 = P(2,1)/(P(1,2)+P(2,1));

% kth moment of sqrt(gamma_cascaded)
mu_cascaded = @(k) kmoment_gamma_cascaded(N, k, K_UR, K_RD);
% kth moment of sqrt(gamma_direct)
mu_direct = @(k) pi1*(bg_cascaded)^(-k/2)/(1-beta)...
    * integral(@(t) t.^(k/2) .* f_gamma_UD(t), tau, Infinity);
% mean(PL_UD*gamma_UD(PL_UD*gamma_UD > tau))/bg_cascaded

mu_gamma_a2g = @(k) kmoment_gamma_a2g(k, mu_cascaded, mu_direct);
%
m_a2g = 4; Q_a2g = 10;
disp('Matching A2G SNR to SSR:...')
[param_a2g, ~, cdf_snr_a2g] = fitdist_SRS(mu_gamma_a2g(1), mu_gamma_a2g(2), m_a2g, Q_a2g);
%% PDF & CDF of A2G SINR
kappa_a2g = param_a2g(1);
omega_a2g = param_a2g(2);
xi_a2g = param_a2g(3);
alp_a2g = param_a2g(5);
Q_a2g = param_a2g(6);

bg_a2g = (1/AWGN)*bg_cascaded;

disp('Calculating A2G SINR...')
cdf_sinr_a2g = @(x) cdf_sinr(Q_a2g, m_a2g, kappa_a2g, omega_a2g, K_ID, bg_ID, K_JD, bg_JD, x/bg_a2g);
cdf_asymp_sinr_a2g = @(x) cdf_asy_sinr(mu_gamma_a2g(1), mu_gamma_a2g(2), m_a2g, Q_a2g, K_ID, bg_ID, K_JD, bg_JD, gamma_ID + gamma_JD + 1, x/bg_a2g);
%% OP - Asymp
PSdB = 10:2.5:50; R_se = 0.5; tau_se = 2^(2*R_se) - 1;

OP_e2e_sim = zeros(1, length(PSdB));
OP_e2e_ana = zeros(1, length(PSdB));
OP_e2e_asymp = zeros(1, length(PSdB));

for iPS = 1:length(PSdB)
    PS = db2pow(PSdB(iPS)); % dBm
    PU = db2pow(PSdB(iPS)); % dBm
    %
    sinr_e2e = min([PS*sinr_g2a; PU*sinr_a2g]);
    
    cdf_sinr_e2e = @(x) 1 - (1-cdf_sinr_g2a(x/PS))*(1-cdf_sinr_a2g(x/PU));

    cdf_asymp_sinr_e2e = @(x) cdf_asymp_sinr_g2a(x/PS) + cdf_asymp_sinr_a2g(x/PU);
    %
    OP_e2e_sim(iPS) = mean(sinr_e2e < tau_se);
    OP_e2e_ana(iPS) = cdf_sinr_e2e(tau_se);
    OP_e2e_asymp(iPS) = cdf_asymp_sinr_e2e(tau_se);
    %
end
%
figure;
semilogy(PSdB, OP_e2e_sim, '-+'); hold on;
semilogy(PSdB, OP_e2e_ana, '-'); hold on;
semilogy(PSdB, OP_e2e_asymp, '--k'); hold on;
axis([min(PSdB) max(PSdB) 1e-12 1]);
%
% filename = sprintf('OP_e2e_sim_N%d_M%d_b%d.mat', N, M, ibeta);
% save(filename, 'OP_e2e_sim');
% filename = sprintf('OP_e2e_ana_N%d_M%d_b%d.mat', N, M, ibeta);
% save(filename, 'OP_e2e_ana');
% filename = sprintf('OP_e2e_asy_N%d_M%d_b%d.mat', N, M, ibeta);
% save(filename, 'OP_e2e_asymp');
%
leg = legend('Simulation', 'Analytical', 'Asymptotic');
leg.ItemTokenSize = [10, 5];
%
xlabel('$P_{\rm S} = P_{\rm U} = P$ [dBm]', 'Interpreter', 'latex');
ylabel('Outage Probability', 'Interpreter', 'latex');