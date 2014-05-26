clear 

left = 1.0;
right = 2.0;

tot_len = right - left;

Nth_order = 6;
Nth_order_init = 3; % We check h-conv from NH_order_init to NH_order
NH_order_init = 3;
NP_order_init = 7; % Order of polynomial basis: 1==classical FEM
Fst = 1; Snd = 2; Trd = 3;

THerr_dom1 = [];  %length/2^Nh_order;
THerr_dom2 = [];  %length/2^Nh_order;
THerr_dom3 = [];  %length/2^Nh_order;
THerr_vec1 = [];
THerr_vec2 = [];
THerr_vec3 = [];

errtol1 = 1e-16;
errtol2 = 1e-16;
errtol3 = 1e-16;

BdyType = 'DN';
deriv_l = 0;
deriv_r = 1;

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = Nth_order_init: Nth_order

  Nth = 2^(ip-1);
  Ndiv = 2^NH_order_init;  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h1 = NP_order_init+Fst;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h1 = [pvec_h1; NP_order_init+Fst];
  end
  
    ueval = zeros(Nth, 2);
    ueval = ANNEvalUserSolution(left, right, Nth, deriv_l, deriv_r, ueval);
    [xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec_h, pvec_h1, Nth, ueval(:, 1), ueval(:, 2));

    [xvec, yvec, zMat_sol] = ANN2DEvalExactSol(left, right, xvec, yvec);

    max_err_th1 = ANN2DCalcErrors(zMat_apx, zMat_sol)

%   [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec_h, pvec_h1, vdl, vdr);
%   [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
%   max_err_h1 = sp1DCalcErrors(yvec_apx, yvec_sol)

  if max_err_th1 > errtol1
    THerr_dom1 = [THerr_dom1;Nth];
    THerr_vec1 = [THerr_vec1; max_err_th1];
  else
    mxerr1 = max_err_th1
    break;
  end
end

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = Nth_order_init: Nth_order

  Nth = 2^(ip-1);
  Ndiv = 2^NH_order_init;  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h2 = NP_order_init+Snd;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h2 = [pvec_h2; NP_order_init+Snd];
  end
  
    ueval = zeros(Nth, 2);
    ueval = ANNEvalUserSolution(left, right, Nth, deriv_l, deriv_r, ueval);
    [xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec_h, pvec_h2, Nth, ueval(:, 1), ueval(:, 2));

    [xvec, yvec, zMat_sol] = ANN2DEvalExactSol(left, right, xvec, yvec);

    max_err_th2 = ANN2DCalcErrors(zMat_apx, zMat_sol)

%   [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec_h, pvec_h1, vdl, vdr);
%   [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
%   max_err_h1 = sp1DCalcErrors(yvec_apx, yvec_sol)

  if max_err_th2 > errtol2
    THerr_dom2 = [THerr_dom2;Nth];
    THerr_vec2 = [THerr_vec2; max_err_th2];
  else
    mxerr3 = max_err_th2
    break;
  end
end

%%%%%%%%%%Loop checking Error for the change of number of elements%%%%%%%%%%%
for ip = Nth_order_init: Nth_order

  Nth = 2^(ip-1);
  Ndiv = 2^NH_order_init;  
  unit_len = tot_len/Ndiv;

  dxvec_h = unit_len;
  pvec_h3 = NP_order_init+Trd;
  
  for j = 2:Ndiv
    dxvec_h = [dxvec_h; unit_len];
    pvec_h3 = [pvec_h3; NP_order_init+Trd];
  end
  
    ueval = zeros(Nth, 2);
    ueval = ANNEvalUserSolution(left, right, Nth, deriv_l, deriv_r, ueval);
    [xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec_h, pvec_h3, Nth, ueval(:, 1), ueval(:, 2));

    [xvec, yvec, zMat_sol] = ANN2DEvalExactSol(left, right, xvec, yvec);

    max_err_th3 = ANN2DCalcErrors(zMat_apx, zMat_sol)

%   [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec_h, pvec_h1, vdl, vdr);
%   [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
%   max_err_h1 = sp1DCalcErrors(yvec_apx, yvec_sol)

  if max_err_th3 > errtol3
    THerr_dom3 = [THerr_dom3;Nth];
    THerr_vec3 = [THerr_vec3; max_err_th3];
  else
    mxerr3 = max_err_th3
    break;
  end
end
% 
% slope1 = ANNEvalLogxLogySlope(THerr_dom1, THerr_vec1)
% slope2 = ANNEvalLogxLogySlope(THerr_dom2, THerr_vec2)
% slope3 = ANNEvalLogxLogySlope(THerr_dom3, THerr_vec3)

figure(2);
  semilogy(THerr_dom1, THerr_vec1, 'bo-','markersize',5, 'linewidth', 2); hold on;
  semilogy(THerr_dom2, THerr_vec2, 'rs--','markersize',5, 'linewidth', 2); hold on;
  semilogy(THerr_dom3, THerr_vec3, 'kv-.','markersize',5, 'linewidth', 2);
    grid on;
    xlabel('Number of modes (|\theta|)');
    ylabel('Discrete L_{\infty} Error');
    pic1 = cat(2, 'Order ', num2str(NP_order_init + Fst));
    pic2 = cat(2, 'Order ', num2str(NP_order_init + Snd));
    pic3 = cat(2, 'Order ', num2str(NP_order_init + Trd));
    
    legend(pic1, pic2, pic3,2);
  hold off;