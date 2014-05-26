left = 0.0;
right = 1.0;

tot_len = right - left;

%%%%%%%%%%Loop checking Error for the change of order of elements%%%%%%%%%%%


NP_order_start = 1;
NP_order_end = 12;

NP_ints = 5;

Perr_dom1 = [];
Perr_vec1 = [];  

for ip = NP_order_start : NP_order_end

  unit_len = tot_len/NP_ints;

  dxvec = unit_len;
  pvec = ip;
  for j = 2:NP_ints
    dxvec = [dxvec;unit_len];
    pvec = [pvec; ip];
  end
  
  Psamples = 5 * ip;
  max_err_p = spDiffStiffnessDD(left, right, dxvec, pvec, 999, 999, Psamples);

  if max_err_p > 5e-15
    Perr_dom1 = [Perr_dom1, ip];
    Perr_vec1 = [Perr_vec1, max_err_p];
  else 
    mxerp05 = max_err_p
    break;
  end
  
end

NP_ints = 10;

Perr_dom2 = [];
Perr_vec2 = [];  

for ip = NP_order_start : NP_order_end

  unit_len = tot_len/NP_ints;

  dxvec = unit_len;
  pvec = ip;
  for j = 2:NP_ints
    dxvec = [dxvec;unit_len];
    pvec = [pvec; ip];
  end
  
  Psamples = 5 * ip;
  max_err_p = spDiffStiffnessDD(left, right, dxvec, pvec, 999, 999, Psamples);

  if max_err_p > 5e-15
    Perr_dom2 = [Perr_dom2, ip];
    Perr_vec2 = [Perr_vec2, max_err_p];
  else 
    mxerp10 = max_err_p  
    break;
  end
  
end

figure(3);

semilogy(Perr_dom1, Perr_vec1, 'bo-','markersize', 6,'linewidth', 3); hold on;
semilogy(Perr_dom2, Perr_vec2, 'rs--','markersize', 6,'linewidth', 3);
    grid on;
    xlabel('order per element');
    ylabel('discrete L-infinite error');
legend('5 elements','10 elements');