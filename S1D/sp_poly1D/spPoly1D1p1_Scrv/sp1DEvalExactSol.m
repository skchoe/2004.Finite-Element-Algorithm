%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);

%   %------------------------------------------------------------------------
%   % f = sin(pi*x); -> u = - pi^(-2) * sin(pi*x), DB(0) = 0, NB(1) = 1/pi
%   fvec = sin(pi*domvec);
%   
%   pi2 = pi*pi;
% 
%   if vdl == 999
%     vLDB = - sin(pi*left) / pi2;
%   else
%     vLDB = vdl;
%   end
%
%   if vnr == 999
%     vRNB = - cos(pi*right) / pi;
%   else
%     vRNB = vnr;
%   end
%   
%   tmp1 = 1/pi * cos(pi*right);
%   tmp2 = 1/(pi2) * sin(pi*left);
%   
%   C = vRNB + tmp1;
%   D = vLDB - left * vRNB + tmp2 - left * tmp1;

%   Ndomvec = size(domvec, 1); 
%   uvec = zeros(Ndomvec, 1);
%   uvec = -1/pi2 * sin(pi * domvec) + C * domvec + D;

%   %------------------------------------------------------------------------
%   % Given K>0, f = K*(K-1)x^(K-2) -> u = x^K, DB(0) = 1, NB(1) = K
%   K = 21;
%   x1vec = power(domvec, K-2);
%   fvec = K * (K-1) * x1vec;
%   
% 
%   if vdl == 999
%     vLDB = left ^ K;
%   else
%     vLDB = vdl;
%   end
% 
%   if vnr == 999
%     vRNB = K * right ^ (K-1);
%   else
%     vRNB = vnr;
%   end
%   
%   tmp1 = right^(K-1);
%   tmp2 = left^K;
%   
%   C = vRNB - K * tmp1;
%   D = vLDB - tmp2 - (vRNB - K * tmp1) * left;

%   Ndomvec = size(domvec, 1);
%   uvec = zeros(Ndomvec, 1);
%   dompwr = power(domvec, K);
%   uvec = dompwr+ C * domvec + D;

%   %------------------------------------------------------------------------
%   % Other function, 2nd derivative, Neumann, Dirichlet BDY's
%   prob_order = 9;
%   u = spPolySCurve(prob_order); % Polynomial wit Moment Basis
% 
%   K = size(u, 1) - 1;
%   fvec = zeros(size(domvec,1), 1);
%   
%   for i = 2:K
%     x1vec = power(domvec, i-2);
%     fvec = fvec + i * (i-1) * u(i+1,1) * x1vec;
%   end
%   
% 
%   vsumL = powersum(left, K, u);
%   if vdl == 999
%     vLDB = vsumL;
%   else
%     vLDB = vdl;
%   end
% 
%   vsumR = powersumDeriv1(right, K, u);
%   if vnr == 999
%     vRNB = vsumR;
%   else
%     vRNB = vnr;
%   end
%   
%   C = vRNB - vsumR;
%   D = vLDB - vsumL - C * left;

%   Ndomvec = size(domvec, 1);
%   uvec = zeros(Ndomvec, 1);
%   uvec = powersum(domvec, K, u) + C * domvec + D;

    %------------------------------------------------------------------------
    % Other function, 2nd derivative, Neumann, Dirichlet BDY's
    prob_order = 9;
    alphas = spPolySCurveLe(prob_order); % Polynomial with Legendre Basis
    Nalpha = size(alphas, 1);
  
    Nxvec = size(xvec, 1);
    uvec = spEvalPolynomialLe(alphas, Nalpha, xvec, Nxvec);
  

    if  BdyType == 'DN'
        vsumL = spEvalPolynomialLe(alphas, Nalpha, left, 1);
        vsumR = spEvalPolynomialLeDeriv1(alphas, Nalpha, right, 1);
        
        C = vBR - vsumR;
        D = vBL - vsumL - C * left;

        yvec_sol = uvec + C * xvec + D;
    
    
    elseif BdyType == 'DD'
        vsumL = spEvalPolynomialLe(alphas, Nalpha, left, 1);
        vsumR = spEvalPolynomialLe(alphas, Nalpha, right, 1);

        d1 = vBR - vBL;
        e1 = vsumR - vsumL;
        dlen = right - left;
    
        C = (d1 - e1) / dlen;
        D = vBR - vsumR - C * right;
    
        yvec_sol = uvec + C * xvec + D;

    else
    end

  %--------------------------------------------------------------  
  figno = 21;
  figure(figno);
      % Analytic solution considering boundary condition
      plot(xvec, yvec_sol, '.-', 'markersize', 5);
      grid on, title('Analytic Solution');    

return

%Subroutines%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outvec = power(iNdomvec, n)

num = size(iNdomvec, 1);
for i=1:num
    outvec(i, 1) = iNdomvec(i, 1) ^ n;
end

return


% n = order == size(u_hat)-1
function outvec = powersum(iNdomvec, n, u_hat)

  outvec = zeros(size(iNdomvec, 1), 1);

  for i = 1:n+1
    outvec = outvec + u_hat(i, 1) * power(iNdomvec, i-1);
    
  end

return

function outvec = powersumDeriv1(iNdomvec, n, u_hat)

  outvec = zeros(size(iNdomvec, 1), 1);

  for i = 2:n+1
    outvec = outvec + (i-1) * u_hat(i, 1) * power(iNdomvec, i-2);
    
  end

return
