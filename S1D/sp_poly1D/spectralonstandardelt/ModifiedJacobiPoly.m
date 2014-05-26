%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ModifiedJacobiPoly.m 
% Evaluation of Modified Jacobi Polynomial
%
%                           (1 - xi)/2    if p == 0
% pi[p](xi) -> psi[p](xi) = (1 + xi)/2    if p == P
%                           (1 - xi)/2 * (1 + xi)/2 * P(1,1)[p-1](xi), if 2 <= p <= P                        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = ModifiedJacobiPoly(degree, x, alpha, beta)

  if degree == 0
    y = 0.5*(1-x);
  elseif degree == 1
    y = 0.5*(1+x);
  else
    y = 0.25 * (1-x).*(1+x) .* JacobiPoly(degree-2, x, alpha, beta);
  end
  
return
