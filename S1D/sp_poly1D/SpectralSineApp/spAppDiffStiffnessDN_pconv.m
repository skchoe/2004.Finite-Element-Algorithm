left = 0.0;
right = 1.0;

tot_len = right - left;

%%%%%%%%%%Loop checking Error for the change of order of elements%%%%%%%%%%%


NP_order_start = 1;
NP_order_end = 15;

NP_ints = 5;

Perr_dom = [];
Perr_vec = [];  

for ip = NP_order_start : NP_order_end

  unit_len = tot_len/NP_ints;

  dxvec = unit_len;
  pvec = ip;
  for j = 2:NP_ints
    dxvec = [dxvec;unit_len];
    pvec = [pvec; ip];
  end
  
  Psamples = 5 * ip;
  max_err_p = spDiffStiffnessDN(left, right, dxvec, pvec, 999, 0, Psamples);

  Perr_dom = [Perr_dom, ip];
  Perr_vec = [Perr_vec, max_err_p];
  
end

figure(3);

  semilogy(Perr_dom, Perr_vec, 'o-','markersize',6, 'linewidth', 3);
    grid on;
    xlabel('order per element');
    ylabel('discrete L-infinite error');

    minerror = min(Perr_vec)