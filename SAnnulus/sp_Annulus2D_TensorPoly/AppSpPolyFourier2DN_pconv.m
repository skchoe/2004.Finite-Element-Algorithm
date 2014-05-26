left = 0.0;
right = 1.0;

tot_len = right - left;

%%%%%%%%%%Loop checking Error for the change of order of elements%%%%%%%%%%%


NP_order_start = 1;
NP_order_end = 10;

NP_ints = 5;

Perr_dom1 = [];
Perr_vec1 = [];  

BdyType = 'DN';
errbound = 1e-15;
errbound2 = 1e-15;

% Boundary vectors preprocessing
deriv_l = 0;    % degree of derivative for left  bdy condition
deriv_r = 1;    % degree of derivative for right bdy condition


Nth = 8;


for ip = NP_order_start : NP_order_end

  unit_len = tot_len/NP_ints;

  dxvec = unit_len;
  pvec = ip;
  for j = 2:NP_ints
    dxvec = [dxvec;unit_len];
    pvec = [pvec; ip];
  end
  
ueval = zeros(Nth, 2);
ueval = ANNEvalUserSolution(left, right, Nth, deriv_l, deriv_r, ueval);

tendeg = 1; % \sigma(r) = r^tendeg; -  should be checked in each routine
sigleft = ANNEvalUserTensor(tendeg, left);
sigright = ANNEvalUserTensor(tendeg, right);

[xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec, pvec, Nth, ueval(:, 1), ueval(:, 2), sigleft, sigright);

[xvec, yvec, zMat_sol] = ANN2DEvalExactSol(left, right, xvec, yvec);

max_err_p1 = ANN2DCalcErrors(zMat_apx, zMat_sol)

%   [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec, pvec, vdl, vdr);
%   [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
%   max_err_p1 = sp1DCalcErrors(yvec_apx, yvec_sol)
% 
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
  
ueval = zeros(Nth, 2);
ueval = ANNEvalUserSolution(left, right, Nth, deriv_l, deriv_r, ueval);

tendeg = 1; % \sigma(r) = r^tendeg; -  should be checked in each routine
sigleft = ANNEvalUserTensor(tendeg, left);
sigright = ANNEvalUserTensor(tendeg, right);

[xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec, pvec, Nth, ueval(:, 1), ueval(:, 2), sigleft, sigright);

[xvec, yvec, zMat_sol] = ANN2DEvalExactSol(left, right, xvec, yvec);

max_err_p2 = ANN2DCalcErrors(zMat_apx, zMat_sol)

%   [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec, pvec, vdl, vdr);
%   [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
%   max_err_p2 = sp1DCalcErrors(yvec_apx, yvec_sol);
  
  if max_err_p2 > errbound2
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
legend('5 Elements','10 Elements');