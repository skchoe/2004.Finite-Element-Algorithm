function  ynodes = spEvalMultiSpan1DNodes(left, right, xnodes, alpha, beta, ug_hat, pvec, map)

  N = size(xnodes, 1);
  ynodes = zeros(N, 1);

  for i = 1:N
    
    x = xnodes(i, 1);

    if i == 1 | i == N 
    
      ynodes(i, 1) = -1 * 1/pi^2 * sin(pi*x); 
      
    else

      eleft  = xnodes(i-1, 1);
      eright = xnodes(i, 1);
      elength = eright - eleft;
      eth    = i-1;

      eorder = pvec(eth, 1);
      gidx = map(eth, :)';

      eug_hat = zeros(1, eorder+1);
      exival = zeros(eorder+1, 1);

      for j = 1:eorder + 1
        eug_hat(1, j) = ug_hat(gidx(j, 1), 1);

        xi = (2*x - (eright + eleft))/elength;   % xi = reparametrization of x to [-1, 1]
        exival(j, 1) = ModifiedJacobiPoly(j-1, xi, alpha, beta);  % Starts from degree 0

      end
      ynodes(i, 1) = eug_hat * exival;
    end

  end


return
