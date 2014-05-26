left = 0.0;
right = 1.0;

tot_len = right - left;

%%%%%%%%%%Loop checking Error for the change of order of elements%%%%%%%%%%%


NP_order_start = 1;
NP_order_end = 9;

NP_ints = 5;

Perr_dom1 = [];
Perr_vec1 = [];  

BdyType = 'DD';
errbound = 5e-15
vdl = 999;
vdr = 999;

for ip = NP_order_start : NP_order_end

  unit_len = tot_len/NP_ints;

  dxvec = unit_len;
  pvec = ip;
  for j = 2:NP_ints
    dxvec = [dxvec;unit_len];
    pvec = [pvec; ip];
  end
  
  [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec, pvec, vdl, vdr);
  [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
  max_err_p1 = sp1DCalcErrors(yvec_apx, yvec_sol)

  if max_err_p1 > errbound
    Perr_dom1 = [Perr_dom1, ip];
    Perr_vec1 = [Perr_vec1, max_err_p1];
  else 
    mxerp05 = max_err_p1
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
  
  [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec, pvec, vdl, vdr);
  [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
  max_err_p2 = sp1DCalcErrors(yvec_apx, yvec_sol)

  if max_err_p2 > errbound
    Perr_dom2 = [Perr_dom2, ip];
    Perr_vec2 = [Perr_vec2, max_err_p2];
  else 
    mxerp10 = max_err_p2  
    break;
  end
  
end

max_err_p1
max_err_p2
figure(3);

semilogy(Perr_dom1, Perr_vec1, 'bo-','markersize', 5,'linewidth', 2); hold on;
semilogy(Perr_dom2, Perr_vec2, 'rs--','markersize', 5,'linewidth', 2);
    grid on;
    xlabel('Order Per Element');
    ylabel('Discrete L_{\infty} Error');
legend('5 elements','10 elements');