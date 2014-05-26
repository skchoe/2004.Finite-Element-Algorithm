%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function crv_img = spEvalUserFunctions(alphas, Nalpha, crv_dom, Nsamp)
% alphas: array of coefficient from order 0 to Nalpha
% crv_dom: vector of domain
% Nsamp:size of crv_dom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uvec, fvec, vLDB, vRNB] = spEvalUserFunctions(left, right, domvec, Nvec, vdl, vnr)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % f = sin(pi*x); -> u = - pi^(-2) * sin(pi*x), DB(0) = 0, NB(1) = 1/pi
  fvec = sin(pi*domvec);
  
  pi2 = pi*pi;

  if vdl == 999
    vLDB = - sin(pi*left) / pi2;
  else
    vLDB = vdl;
  end

  if vnr == 999
    vRNB = - cos(pi*right) / pi;
  else
    vRNB = vnr;
  end
  
  tmp1 = 1/pi * cos(pi*right);
  tmp2 = 1/(pi2) * sin(pi*left);
  
  C = vRNB + tmp1;
  D = vLDB - left * vRNB + tmp2 - left * tmp1;

  uvec = zeros(Nvec, 1);
  uvec = -1/pi2 * sin(pi * domvec) + C * domvec + D;

%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Given K>0, f = K*(K-1)x^(K-2) -> u = x^K +1, DB(0) = 1, NB(1) = K
%   K = 3;
%   x1vec = domvec^(K-2);
%   fvec = K * (K-1) * x1vec;
%   
%   x2vec = x1vec^2;
%   uvec = x2vec + 1;
%   
%   vLDB = left^K + 1;
%   vRNB = K * right^(K-1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Other function, 2nd derivative, Neumann, Dirichlet BDY's
  

  
return
