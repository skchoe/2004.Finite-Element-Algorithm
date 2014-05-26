%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [uvec, fvec, vLDB, vRDB] = spEvalUserFunctionsDD(left, right, domvec, vdl, vdr)
% alphas: array of coefficient from order 0 to Nalpha
% crv_dom: vector of domain
% Nsamp:size of crv_dom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fvec, vLDB, vRDB] = spEvalUserFunctionsDD(left, right, domvec, vdl, vdr)

%   %------------------------------------------------------------------------------------------------------------------------------
  prob_order = 9;
  alphas = spPolySCurveLe(prob_order); % Polynomial with Legendre Basis
  Nalpha = size(alphas, 1);
 
  Ndomvec = size(domvec, 1);
  fvec = spEvalPolynomialLeDeriv2(alphas, Nalpha, domvec, Ndomvec);
  
  vsumL = spEvalPolynomialLe(alphas, Nalpha, left, 1);
  if vdl == 999
    vLDB = vsumL;
  else
    vLDB = vdl;
  end

  vsumR = spEvalPolynomialLe(alphas, Nalpha, right, 1);
  if vdr == 999
    vRDB = vsumR;
  else
    vRDB = vdr;
  end
  
%   domlen = right - left;
%   ulen = tmpr - tmpl;
%   bdlen = vRDB - vLDB;
%   
%   C = (bdlen - ulen) / domlen;
%   D = vLDB - tmpl - C*left;
% 
%   Ndomvec = size(domvec, 1);
%   uvec = zeros(Ndomvec, 1);
%   uvec = -1/pi2 * sin(pi * domvec) + C * domvec + D;
% 
  
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