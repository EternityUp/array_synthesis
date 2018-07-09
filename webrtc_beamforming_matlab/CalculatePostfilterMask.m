function mask = CalculatePostfilterMask...
    (interf_cov_mat,rpsiw,ratio_rxiw_rxim,rmw_r,eig_m,kCutOffConstant)
%¼ÆËãºóÖÃÂË²¨Æ÷ÑÚ±ÎÖµ
rpsim = conj(eig_m) * interf_cov_mat * eig_m.';
if ( rpsim > 0.0 )
    ratio = rpsiw / rpsim;
end
numerator = 1.0 - kCutOffConstant;
if ( rmw_r > 0.0 )
    numerator = 1.0 - min(kCutOffConstant, ratio / rmw_r);
end
denominator = 1.0 - kCutOffConstant;
if ( ratio_rxiw_rxim > 0.0 )
    denominator = 1.0 - min(kCutOffConstant, ratio / ratio_rxiw_rxim);
end
mask = real( numerator / denominator );