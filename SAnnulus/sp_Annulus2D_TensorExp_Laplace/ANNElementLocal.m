
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANNLocalMatsPolyFourier2D.m Local Matrix of Spectral Element Method in 1D
%
% ANNLocalMatsPolyFourier2D(P, alpha, beta, w, z, Jn)
%   P: Order of approximation
%   alpha, beta: of Modified Jacobi Polynomial
%   w, z: vectors from Gaussian Quadrature
%   Jn: Jacobian of mapping from [-1,1] to [a, b]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M1, M2, M3] = ANNElementLocal(Pn, alpha, beta, z, w, xiz, Jn)

    % Computing L and fvector
    
    for i = 0:Pn
        Lv(:, i+1) =  ModifiedJacobiPoly(i, Pn, z, alpha, beta);
        Lpv(:, i+1) = ModifiedJacobiPolyDerivative(i, Pn, z, alpha, beta);
    end
 
    tiz = ANNEvalUserTensor(xiz);
    
    for p = 0:Pn

        for j = 0:Pn
      
            M1(p+1, j+1) = w' * (xiz .*  tiz .* Lv(:, p+1)  .* Lpv(:, j+1));
            
            M2(p+1, j+1) = w' * (xiz .* xiz .*  tiz .* Lpv(:, p+1) .* Lpv(:, j+1)) / Jn;
            
            M3(p+1, j+1) = w' * (tiz .* Lv(:, p+1) .* Lv(:, j+1)) * Jn;
            
            if abs(M1(p+1, j+1)) < 1.0e-15, M1(p+1, j+1) = 0.0; end
            if abs(M2(p+1, j+1)) < 1.0e-15, M2(p+1, j+1) = 0.0; end
            if abs(M3(p+1, j+1)) < 1.0e-15, M3(p+1, j+1) = 0.0; end
       end
  
    end

return