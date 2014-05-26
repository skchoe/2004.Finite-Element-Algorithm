left = 0.0;
right = 1.0;

tot_len = right - left;

NH_start = 2;
NH_order = 15;
NH_order_init = 1; % Order of polynomial basis: 1==classical FEM

Herr_dom1 = [];  %length/2^Nh_order;
Herr_dom2 = [];  %length/2^Nh_order;
Herr_dom3 = [];  %length/2^Nh_order;
Herr_vec1 = [];
Herr_vec2 = [];
Herr_vec3 = [];


%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = NH_start : NH_order

  Ndiv = 2^(ip-1);  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h1 = NH_order_init;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h1 = [pvec_h1; NH_order_init];
  end

  % Number of samples to draw graph or find error
  Hsamples = Ndiv + 1;
  
  max_err_h1 = spDiffStiffnessDN(left, right, dxvec_h, pvec_h1, 999, 0, Hsamples);

  if max_err_h1 > 5e-15
    Herr_dom1 = [Herr_dom1; unit_len];
    Herr_vec1 = [Herr_vec1; max_err_h1];
  else
    mxerr1 = max_err_h1
    break;
  end
end

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = 1 : NH_order

  Ndiv = 2^(ip-1);  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h2 = NH_order_init + 2;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h2 = [pvec_h2; NH_order_init + 2];
  end

  % Number of samples to draw graph or find error
  Hsamples = Ndiv + 1;
  
  max_err_h2 = spDiffStiffnessDN(left, right, dxvec_h, pvec_h2, 999, 0, Hsamples);

  if max_err_h2 > 5e-15
    Herr_dom2 = [Herr_dom2; unit_len];
    Herr_vec2 = [Herr_vec2; max_err_h2];  
  else
    mxerr2 = max_err_h2
    break;
  end
end

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = 1 : NH_order

  Ndiv = 2^(ip-1);  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h3 = NH_order_init + 5;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h3 = [pvec_h3; NH_order_init + 5];
  end

  % Number of samples to draw graph or find error
  Hsamples = Ndiv + 1;
  
  max_err_h3 = spDiffStiffnessDN(left, right, dxvec_h, pvec_h3, 999, 0, Hsamples)

  if max_err_h3 > 5e-15
    Herr_dom3 = [Herr_dom3; unit_len];
    Herr_vec3 = [Herr_vec3; max_err_h3];  
  else
    mxerr3 = max_err_h3
    break;
  end
  
end

slope1 = spEvalLogxlogySlope(Herr_dom1, Herr_vec1);
slope2 = spEvalLogxlogySlope(Herr_dom2, Herr_vec2);
slope3 = spEvalLogxlogySlope(Herr_dom3, Herr_vec3);

figure(2);
  loglog(Herr_dom1, Herr_vec1, 'bo-','markersize',6, 'linewidth', 3); hold on;
  loglog(Herr_dom2, Herr_vec2, 'rs--','markersize',6, 'linewidth', 3); hold on;
  loglog(Herr_dom3, Herr_vec3, 'kv-.','markersize',6, 'linewidth', 3);
    grid on;
    xlabel('element size(h)');
    ylabel('discrete L-infinite error');
    pic1 = cat(2, 'order 1 of slope ', num2str(slope1));
    pic2 = cat(2, 'order 3 of slope ', num2str(slope2));
    pic3 = cat(2, 'order 6 of slope ', num2str(slope3));
    
    legend(pic1, pic2, pic3,2);
  hold off;