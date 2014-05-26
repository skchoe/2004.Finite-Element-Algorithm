left = 0.0;
right = 1.0;

tot_len = right - left;

NH_order = 3;
NH_order_init = 1; % Order of polynomial basis: 1==classical FEM

Herr_dom = [];  %length/2^Nh_order;
Herr_vec1 = [];
Herr_vec3 = [];
Herr_vec5 = [];


%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = 1 : NH_order

  Ndiv = 2^(NH_order - ip);  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h1 = NH_order_init;
  pvec_h3 = NH_order_init + 3;
  pvec_h5 = NH_order_init + 5;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h1 = [pvec_h1; NH_order_init];
    pvec_h3 = [pvec_h3; NH_order_init + 3];
    pvec_h5 = [pvec_h5; NH_order_init + 5];
  end

  % Number of samples to draw graph or find error
  Hsamples = Ndiv + 1;
  
  max_err_h1 = spDiffStiffnessDD(left, right, dxvec_h, pvec_h1, 999, 999, Hsamples);
  max_err_h3 = spDiffStiffnessDD(left, right, dxvec_h, pvec_h3, 999, 999, Hsamples);
  max_err_h5 = spDiffStiffnessDD(left, right, dxvec_h, pvec_h5, 999, 999, Hsamples);

  Herr_dom = [Herr_dom, unit_len];
  Herr_vec1 = [Herr_vec1, max_err_h1];
  Herr_vec3 = [Herr_vec3, max_err_h3];  
  Herr_vec5 = [Herr_vec5, max_err_h5];  
end

figure(2);

  loglog(Herr_dom, Herr_vec1, 'o-','markersize',6, 'linewidth', 3); hold on;
  loglog(Herr_dom, Herr_vec3, 's:','markersize',6, 'linewidth', 3); hold on;
  loglog(Herr_dom, Herr_vec5, 'v-.','markersize',6, 'linewidth', 3);
    grid on;
    xlabel('element size(h)');
    ylabel('discrete L-infinite error');
  hold off;

