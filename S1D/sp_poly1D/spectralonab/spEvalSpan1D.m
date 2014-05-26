%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spEvalSpan.m Evaluation of vector input of Linear combination of basis in 1D
%
% spEvalSpan(P, alpha, beta, left, right, noelts)
%   P: Order of approximation
%   alpha, beta: of Modified Jacobi Polynomial
%   left, right: Interval of domain
%   noelts: number of element interval in domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xvec, yvec] = spEvalSpan1D(P, alpha, beta, left, right, noelts, u_hat)

  h = (right - left)/noelts;

  xvec = [left:h:right]';
  sx = 2/(right-left) * (xvec - 0.5*(left + right));
  
  MjcbPoly = ModifiedJacobiPoly(0, sx, alpha, beta);  % Starts from degree 0
  if P > 0
    for l = 1:P   %Degree of Jacobi polynomial in sxi
      MjcbPolyDeg = ModifiedJacobiPoly(l, sx, alpha, beta);
      MjcbPoly = [MjcbPoly, MjcbPolyDeg];
    end
  end
  
  yvec = MjcbPoly * u_hat;
  
return