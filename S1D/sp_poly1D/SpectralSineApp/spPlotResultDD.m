%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function spPlotResultDD(intervals, left, right, alpha, beta, dxvec, pvec, ug_hat, map, figno)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function max_error = spPlotResultDD(intervals, left, right, alpha, beta, dxvec, pvec, Lg, ug_hat, map, figno)

  N_hat = size(ug_hat, 1);
  
  length = right - left;
  h = length/intervals;

  % vector for domain in plot routine
  xvec = [left:h:right]';
  
  % array of nodes
  N_node = size(dxvec,1);
  xnodes = zeros(N_node, 1);
  xnodes(1, 1) = left;
  for i = 1:N_node
    xnodes(i+1, 1) = xnodes(i, 1) + dxvec(i, 1);
  end

  
  ynodes = spEvalMultiSpan1DNodes(left, right, xnodes, alpha, beta, ug_hat, pvec, map);
  % output vector for input values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %xvec, eleftvec, erightvec, alpha, beta, ug_hat, pvec, map, ethvec
  yvec = spEvalMultiSpan1D(left, right, dxvec, xvec, alpha, beta, ug_hat, pvec, map);
  figure(figno);
  
  subplot(3, 1, 1);
    % solution curve
    plot(xvec, yvec, '.-','markersize',13);
    title('Numerical Solution');
    hold on

    % elements and orders (elements, orders ploting)
    N_elt = size(dxvec, 1);
    tmp = left; xpos = left;
    for i = 1:N_elt
      x1 = tmp + dxvec(i, 1);
      xpos = [xpos, x1];
      tmp = x1;
    end
    ypos = zeros(N_elt + 1, 1);
    plot(xpos, ypos, 'rx-','markersize',10);
    grid on, title('Approximation');
    legend('Solution','Elements with Orders', 0);

    % text plotting showing <order> of each element
    for j = 1:N_elt
      tx = 0.5 * (xpos(1, j)+xpos(1, j+1)) + 0.01*length;
      ty = -0.05*length;
      todr = num2str(pvec(j,1));
      text(tx, ty, todr);
    end
%
    subplot(3, 1, 2);
    % plot of analytic solution
    yv = zeros(intervals+1, 1);
  
    yv = -1/(pi^2) * sin(pi * xvec);

    max_error = max(abs(yvec - yv));

    plot(xvec, yv, '.-', 'markersize',13);
    grid on, title('Analytic Solution');    
    hold on
    plot(xpos, ypos, 'rx-', 'markersize',10);
    
  % plotting of global stiffness matrix Lg
    subplot(3, 1, 3);
    zerorange = 1e-13;
    spSetSpy(Lg, zerorange, figno, 1);


return