
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function yvec = spEvalMultiSpan1D(xvec, eleftvec, erightvec, alpha, beta, ug_hat, pvec, map, ethvec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function yvec = spEvalMultiSpan1D(xvec, eleftvec, erightvec, alpha, beta, ug_hat, pvec, map, ethvec)

function yvec = spEvalMultiSpan1D(left, right, dxvec, xvec, alpha, beta, ug_hat, pvec, map)
% element range(left, right), order on the element for each samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  [eleftvec, erightvec, ethvec] = spFindIntervalByValue(left, right, dxvec, xvec);

  N = size(xvec, 1);
  yvec = zeros(N, 1);

  for i = 1:N
    x = xvec(i, 1);
    eleft  = eleftvec(i, 1);
    eright = erightvec(i, 1);
    elength = eright - eleft;
    eth    = ethvec(i, 1);

    eorder = pvec(eth, 1);
    gidx = map(eth, :)';

    eug_hat = zeros(1, eorder+1);
    exival = zeros(eorder+1, 1);

    for j = 1:eorder + 1
      eug_hat(1, j) = ug_hat(gidx(j, 1), 1);

      xi = (2*x - (eright + eleft))/elength;   % xi = reparametrization of x to [-1, 1]
      exival(j, 1) = ModifiedJacobiPoly(j-1, xi, alpha, beta);  % Starts from degree 0

%       j_order = j
%       xi_xvalnorm = xi
%       basiseval = exival(j, 1)
%       
    end
    yvec(i, 1) = eug_hat * exival;

  end
%   xvec
%  
%   eug_hat

return

