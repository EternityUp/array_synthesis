function eeratio = ComputeLogEnergyEntropyRatio(pf,ini_ind,end_ind,alpha_log)

pf_seg = pf(ini_ind:end_ind);
sum_pf_seg = sum(pf_seg);

pf_pdf = pf_seg / sum_pf_seg;
pf_pdf = max(pf_pdf,1e-6);

entropy = sum(-pf_pdf .* log2(pf_pdf));

log_spectrum = log10( 1 + sum_pf_seg / alpha_log );

eeratio = log_spectrum / entropy;




