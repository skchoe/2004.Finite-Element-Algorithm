
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spLocalMats1D.m Local Matrix of Spectral Element Method in 1D
%
% spLocalMats1D(P, alpha, beta, w, z, fzs, Jn)
%   P: Order of approximation
%   alpha, beta: of Modified Jacobi Polynomial
%   w, z: vectors from Gaussian Quadrature
%   fzs: vector of evaluation of RHS force
%   Jn: Jacobian of mapping from [-1,1] to [a, b]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, fvector] = spLocalMats1D(P, alpha, beta, z, w, kFun, fFun, xiz, Jn)

  tiz = feval(kFun, xiz);
  fiz = feval(fFun, xiz);
  
  for p = 1:P+1
      
    for q = 1:P+1

      mjcbvec1 = ModifiedJacobiPolyDerivative(p-1, P, z, alpha, beta);
      mjcbvec2 = ModifiedJacobiPolyDerivative(q-1, P, z, alpha, beta);
      % sum up for innerproduct w.uz - LHS

      % Note that each element has negative sign/
      L(p, q) = - w' * (tiz .* mjcbvec1 .* mjcbvec2) / Jn;

    end

    jv = ModifiedJacobiPoly(p-1, P, z, alpha, beta);
    fz = fiz .* jv;

    fvector(p, 1) = Jn * w' * fz;
  
  end
  
return