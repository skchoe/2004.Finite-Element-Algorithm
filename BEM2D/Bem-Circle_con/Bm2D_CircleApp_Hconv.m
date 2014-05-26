R = 1.0;

Np = 4;    % Order of quadrature

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NH_order = 10;
NH_order_init = 5;

Herr_dom1 = [];  %length/2^Nh_order;
Herr_vec1 = [];

errtol = 1e-15;

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = NH_order_init: NH_order

  Ne = 2^(ip-1);
  unit_len = 2*pi*R/Ne;

  dxvec_h = unit_len;
  
  for j = 2:Ne
    dxvec_h = [dxvec_h; unit_len];
end
 
  us = zeros(Ne, 1); % Left circle boundary values
  qs = zeros(Ne, 1); % Left circle normal values

  usol = zeros(Ne, 1);    % Analytic solution
  [usol, qs] = Bm2D_CircleExact(Ne, R);
  us = Bm2D_Circle(Np, R, Ne, us, qs);
  max_err_h1 = norm(usol - us)

  if max_err_h1 > errtol
    Herr_dom1 = [unit_len; Herr_dom1];
    Herr_vec1 = [max_err_h1;Herr_vec1];
  else
    mxerr1 = max_err_h1
    break;
  end
end


slope1 = Bm2D_EvalSlope(Herr_dom1, Herr_vec1)

figure(2);
  loglog(Herr_dom1, Herr_vec1, 'kv-.','markersize',5, 'linewidth', 2);
    grid on;
    xlabel('Element Size (h)');
    ylabel('Discrete L_{\infty} Error');
%     pic3 = cat(2, 'Order ', num2str(NP_order_init));
    
%     legend(pic1, pic2, pic3,2);
  hold off;