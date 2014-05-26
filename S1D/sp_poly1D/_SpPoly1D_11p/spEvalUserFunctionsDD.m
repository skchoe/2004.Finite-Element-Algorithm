%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [uvec, fvec, vLDB, vRDB] = spEvalUserFunctionsDD(left, right, domvec, Nvec, vdl, vdr)
% alphas: array of coefficient from order 0 to Nalpha
% crv_dom: vector of domain
% Nsamp:size of crv_dom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uvec, fvec, vLDB, vRDB] = spEvalUserFunctionsDD(left, right, domvec, Ndomvec, vdl, vdr)

%   %------------------------------------------------------------------------------------------------------------------------------
%   % f = sin(pi*x); -> u = - pi^(-2) * sin(pi*x), DB(0) = 0, NB(1) = 1/pi
%   fvec = sin(pi*domvec);
%   
%   pi2 = pi*pi;
% 
%   tmpl = - 1/pi2 * sin(pi*left);
%   if vdl == 999
%     vLDB = tmpl;
%   else
%     vLDB = vdl;
%   end
% 
%   tmpr = - 1/pi2 * sin(pi*right);
%   if vdr == 999
%     vRDB = tmpr;
%   else
%     vRDB = vdr;
%   end
%   
%   domlen = right - left;
%   ulen = tmpr - tmpl;
%   bdlen = vRDB - vLDB;
%   
%   C = (bdlen - ulen) / domlen;
%   D = vLDB - tmpl - C*left;
% 
%   uvec = zeros(Ndomvec, 1);
%   uvec = -1/pi2 * sin(pi * domvec) + C * domvec + D;

  
%   %------------------------------------------------------------------------------------------------------------------------------
%   % Given K>0, f = K*(K-1)x^(K-2) -> u = x^K, DB(0) = 1, NB(1) = K
%   K = 21;
%   x1vec = power(domvec, K-2);
%   fvec = K * (K-1) * x1vec;
%   
%   tmpl = left^K; 
%   if vdl == 999
%     vLDB = tmpl;
%   else
%     vLDB = vdl;
%   end
% 
%   tmpr = right^K; 
%   if vdr == 999
%     vRDB = tmpr;
%   else
%     vRDB = vdr;
%   end
%   
%   domlen = right - left;
%   ulen = tmpr - tmpl;
%   bdlen = vRDB - vLDB;
%   
%   C = (bdlen - ulen) / domlen;
%   D = vLDB - tmpl - C*left;
% 
%   uvec = zeros(Ndomvec, 1);
%   uvec = power(domvec, K) + C * domvec + D;
%   

%   %------------------------------------------------------------------------------------------------------------------------------
%   % Given K>0, f = K*(K-1)x^(K-2) -> u = x^K, DB(0) = 1, NB(1) = K
%   u = 1.0e+002 * [...
%                   0;...
%                   0;...
%                   0;...
%                   0;...
%                   0;...
%    1.26000000000018;...
%   -4.20000000000046;...
%    5.40000000000049;...
%   -3.15000000000024;...
%    0.70000000000005];  %U of order =9
% %   u = [           0;... 
% %                   0;...
% %                   0;...
% %                   0;...
% %   34.99999999999997;...
% %  -83.99999999999982;...
% %   69.99999999999980;...
% %  -19.99999999999993];     %Order of U=7
% 
% %   prob_order = 9;
% %   u = spPolySCurve(prob_order)
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
%   vsumR = powersum(right, K, u);
%   if vdr == 999
%     vRDB = vsumR;
%   else
%     vRDB = vdr;
%   end
%   
%   d1 = vRDB - vLDB;
%   e1 = vsumR - vsumL;
%   dlen = right - left;
%   
%   C = (d1 - e1) / dlen;
%   D = vRDB - vsumR - C * right;
% 
%   uvec = zeros(Ndomvec, 1);
%   dompwr = power(domvec, K);
%   uvec = powersum(domvec, K, u) + C * domvec + D;
%   
  %------------------------------------------------------------------------------------------------------------------------------
  % Other function(S curve based on Legendre basis, 2nd derivative, Neumann, Dirichlet BDY's
  prob_order = 9;
  alphas = spPolySCurveLe(prob_order); % Polynomial with Legendre Basis
  Nalpha = size(alphas, 1);
  
  fvec = spEvalPolynomialLeDeriv2(alphas, Nalpha, domvec, Ndomvec);
  uvec_tmp = spEvalPolynomialLe(alphas, Nalpha, domvec, Ndomvec);
  
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
  
  d1 = vRDB - vLDB;
  e1 = vsumR - vsumL;
  dlen = right - left;
  
  C = (d1 - e1) / dlen;
  D = vRDB - vsumR - C * right;

  uvec = uvec_tmp + C * domvec + D;
  
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



%   u = [           0;... 
%                   0;...
%                   0;...
%                   0;...
%   34.99999999999997;...
%  -83.99999999999982;...
%   69.99999999999980;...
%  -19.99999999999993];     %Order of U=7
%   u = 1.0e+002 * [...
%                   0;...
%                   0;...
%                   0;...
%   -0.00000000000004;...
%                   0;...
%                   0;...
%    2.09999999999450;...
%   -7.19999999998155;...
%    9.44999999997615;...
%   -5.59999999998603;...
%    1.25999999999688]        %Order of U=10, double precision.
%   u = 1.0e+005 * [...
%                   0;...
%                   0;...
%                   0;...
%    0.00000000000001;...
%                   0;...
%    0.00000000000009;...
%                   0;...
%                   0;...
%    0.06434999923923;...
%   -0.40039999531433;...
%    1.08107998745801;...
%   -1.63799998113657;...
%    1.50149998281703;...
%   -0.83159999053469;...
%    0.25739999708399;...
%   -0.03431999961279];     %Order of U=15


% u = 1.0e+007 * [...
%                   0;...
%                   0;...
%                   0;...
%   -0.00000000251948;...
%                   0;...
%    0.00000000021671;...
%                   0;...
%    0.00000000000874;...
%                   0;...
%   -0.00000000001827;...
%                   0;...
%    0.03530657231169;...
%   -0.32359796531411;...
%    1.34402251283651;...
%   -3.32774081961163;...
%    5.43487703729389;...
%   -6.11382090044690;...
%    4.79487199011675;...
%   -2.58757605933920;...
%    0.91922896291630;...
%   -0.19405167939687;...
%    0.01848045148534];         %Order of U=21
