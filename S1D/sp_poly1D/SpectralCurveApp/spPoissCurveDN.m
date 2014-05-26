%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function spPoissCurveDN(intervals, left, right, alpha, beta, dxvec, pvec, bdleft, bdright, ug_hat, map, figno)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function max_error = spPoissCurveDN(intervals, left, right, alpha, beta, dxvec, pvec, bdleft, bdright, Lg, ug_hat, map, figno, exactsolvec)

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

  
  ynodes = spEvalMultiSpan1DNodes(left, right, xnodes,alpha, beta, ug_hat, pvec, map);
  % output vector for input values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %xvec, eleftvec, erightvec, alpha, beta, ug_hat, pvec, map, ethvec
  yvec = spEvalMultiSpan1D(left, right, dxvec, xvec, alpha, beta, ug_hat, pvec, map);
  figure(figno);
  
%   subplot(2, 1, 1);
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

  diffvec = abs(yvec - exactsolvec);
  max_error = max(diffvec);
  

  figure(figno+1);
    plot(xvec, diffvec, 'r.-','markersize',10);
    grid on, title('Change of differences');
    ylabel('Difference of approximation');
  
% plotting of global stiffness matrix Lg
%   subplot(2, 1, 2);
%   zerorange = 1e-13;
%   spSetSpy(Lg, zerorange, figno, 1);


return