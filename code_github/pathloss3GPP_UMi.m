function pathloss = pathloss3GPP_UMi(Gr,Gt,fc,x,refDistance)
    if (nargin < nargin('pathloss3GPP_UMi'))
        refDistance = 1;
    end
%     PLdB = - Gr - Gt + 22.7 + 26*log10(fc) - 36.7*log10(x);
%     out = db2pow(PLdB);
    
    antennuationCoeff = 186.2087*(fc)^2.6*10^(-Gr/10-Gt/10);
    pathlosExp = 3.67;
    pathloss = antennuationCoeff * (refDistance./x).^pathlosExp;    
end