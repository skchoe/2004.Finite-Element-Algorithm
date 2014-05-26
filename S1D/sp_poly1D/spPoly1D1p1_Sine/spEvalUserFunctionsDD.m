%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [uvec, fvec, vLDB, vRDB] = spEvalUserFunctionsDD(left, right, domvec, vdl, vdr)
% alphas: array of coefficient from order 0 to Nalpha
% crv_dom: vector of domain
% Nsamp:size of crv_dom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fvec, vLDB, vRDB] = spEvalUserFunctionsDD(left, right, domvec, vdl, vdr)

  %------------------------------------------------------------------------
  % f = sin(pi*x); -> u = - pi^(-2) * sin(pi*x), DB(0) = 0, NB(1) = 1/pi
  fvec = sin(pi*domvec);
  
  pi2 = pi*pi;

  if vdl == 999
    vLDB = - sin(pi*left) / pi2;
  else
    vLDB = vdl;
  end

  if vdr == 999
    vRDB = - sin(pi*right) / pi2;
  else
    vRDB = vdr;
  end
  
%   tmp1 = 1/pi * cos(pi*right);
%   tmp2 = 1/(pi2) * sin(pi*left);
%   
%   C = vRNB + tmp1;
%   D = vLDB - left * vRNB + tmp2 - left * tmp1;
% 
%   Ndomvec = size(domvec, 1); 
%   uvec = zeros(Ndomvec, 1);
%   uvec = -1/pi2 * sin(pi * domvec) + C * domvec + D;
  
  
  
% %   %------------------------------------------------------------------------------------------------------------------------------
%   S-Curve example:
%   prob_order = 9;
%   alphas = spPolySCurveLe(prob_order); % Polynomial with Legendre Basis
%   Nalpha = size(alphas, 1);
%  
%   Ndomvec = size(domvec, 1);
%   fvec = spEvalPolynomialLeDeriv2(alphas, Nalpha, domvec, Ndomvec);
%   
%   vsumL = spEvalPolynomialLe(alphas, Nalpha, left, 1);
%   if vdl == 999
%     vLDB = vsumL;
%   else
%     vLDB = vdl;
%   end
% 
%   vsumR = spEvalPolynomialLe(alphas, Nalpha, right, 1);
%   if vdr == 999
%     vRDB = vsumR;
%   else
%     vRDB = vdr;
%   end
%   
% 
% % 


return



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