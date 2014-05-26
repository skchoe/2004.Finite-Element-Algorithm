left = 0.0;
right = 1.0;

tot_len = right - left;

NH_order = 9;
NH_order_init = 4;
NP_order_init = 2; % Order of polynomial basis: 1==classical FEM
Fst = 1; Snd = 2; Trd = 3;

Herr_dom1 = [];  %length/2^Nh_order;
Herr_dom2 = [];  %length/2^Nh_order;
Herr_dom3 = [];  %length/2^Nh_order;
Herr_vec1 = [];
Herr_vec2 = [];
Herr_vec3 = [];

errtol = 1e-13;
%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = NH_order_init: NH_order

  Ndiv = 2^(ip-1);  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h1 = NP_order_init+Fst;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h1 = [pvec_h1; NP_order_init+Fst];
  end
  
  max_err_h1 = sp1DSolvePoisson('DN', left, right, dxvec_h, pvec_h1, 999, 0);

  if max_err_h1 > errtol
    Herr_dom1 = [unit_len; Herr_dom1];
    Herr_vec1 = [max_err_h1;Herr_vec1];
  else
    mxerr1 = max_err_h1
    break;
  end
end

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = NH_order_init : NH_order

  Ndiv = 2^(ip-1);  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h2 = NP_order_init + Snd;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h2 = [pvec_h2; NP_order_init + Snd];
  end

  max_err_h2 = sp1DSolvePoisson('DN', left, right, dxvec_h, pvec_h2, 999, 0);

  if max_err_h2 > errtol
    Herr_dom2 = [unit_len;Herr_dom2 ];
    Herr_vec2 = [max_err_h2;Herr_vec2];  
  else
    mxerr2 = max_err_h2
    break;
  end
end

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = NH_order_init : NH_order

  Ndiv = 2^(ip-1);
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h3 = NP_order_init + Trd;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h3 = [pvec_h3; NP_order_init + Trd];
  end

  max_err_h3 = sp1DSolvePoisson('DN', left, right, dxvec_h, pvec_h3, 999, 0);

  if max_err_h3 > errtol
    Herr_dom3 = [unit_len;Herr_dom3 ];
    Herr_vec3 = [max_err_h3;Herr_vec3];  
  else
    mxerr3 = max_err_h3
    break;
  end
  
end

max_err_h1
max_err_h2
max_err_h3
slope1 = spEvalLogxLogySlope(Herr_dom1, Herr_vec1)
slope2 = spEvalLogxLogySlope(Herr_dom2, Herr_vec2)
slope3 = spEvalLogxLogySlope(Herr_dom3, Herr_vec3)

figure(2);
  loglog(Herr_dom1, Herr_vec1, 'bo-','markersize',5, 'linewidth', 2); hold on;
  loglog(Herr_dom2, Herr_vec2, 'rs--','markersize',5, 'linewidth', 2); hold on;
  loglog(Herr_dom3, Herr_vec3, 'kv-.','markersize',5, 'linewidth', 2);
    grid on;
    xlabel('Element Size (h)');
    ylabel('Discrete L_{\infty} Error');
    pic1 = cat(2, 'Order ', num2str(NP_order_init + Fst));
    pic2 = cat(2, 'Order ', num2str(NP_order_init + Snd));
    pic3 = cat(2, 'Order ', num2str(NP_order_init + Trd));
    
    legend(pic1, pic2, pic3,2);
  hold off;