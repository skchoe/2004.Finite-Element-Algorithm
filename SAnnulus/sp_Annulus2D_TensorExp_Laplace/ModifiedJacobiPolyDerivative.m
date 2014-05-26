%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ModifiedJacobiPolyDerivative.m 
% Evaluation of Derivative of Modified Jacobi Polynomial
%
%                                       (1 - xi)/2    if p == 0
% pi[p](xi) -> psi[p](xi) = (1 - xi)/2 * (1 + xi)/2 * P(1,1)_{p-1}(xi), if 1<= p <= P-1
%                                       (1 + xi)/2   if p == P
%                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = ModifiedJacobiPolyDerivative(degree, maxdegree, x, alpha, beta)

  if degree == 0
    y = -0.5*ones(size(x));
  elseif degree == maxdegree
    y = 0.5*ones(size(x));
  else
    y = -0.5*x.*JacobiPoly(degree-1, x, alpha, beta) + ...
        0.25*(1-x).*(1+x).*JacobiPolyDerivative(degree-1, x, alpha, beta);
  end
return
