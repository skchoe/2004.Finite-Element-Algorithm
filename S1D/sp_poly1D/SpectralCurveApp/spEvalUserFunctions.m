%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function crv_img = spEvalUserFunctions(alphas, Nalpha, crv_dom, Nsamp)
% alphas: array of coefficient from order 0 to Nalpha
% crv_dom: vector of domain
% Nsamp:size of crv_dom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uvec, fvec, vLDB, vRNB] = spEvalUserFunctions(domvec, Nvec)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % f = sin(pi*x); -> u = - pi^(-2) * sin(pi*x), DB(0) = 0, NB(1) = 1/pi
  fvec = sin(pi*domvec);
  
  pi2 = pi*pi;
  uvec = - fvec /pi2;

  vLDB = 0;
  vRNB = 1/pi;
  
  
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Given K>0, f = K*(K-1)x^(K-2) -> u = x^K +1, DB(0) = 1, NB(1) = K
%   K = 3;
%   x1vec = domvec^(K-2);
%   fvec = K * (K-1) * x1vec;
%   
%   x2vec = x1vec^2;
%   uvec = x2vec + 1;
%   
%   vLDB = 1;
%   vRNB = K;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Other function, 2nd derivative, Neumann, Dirichlet BDY's
  

  
return
