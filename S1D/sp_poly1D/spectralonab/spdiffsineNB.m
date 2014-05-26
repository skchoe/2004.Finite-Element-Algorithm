%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffsineNB.m Spectral Element Solver on 1D
%
%   max_error = spdiffsineNB(left, right, samples, P, vdl, vnr, viz)
%   left, right : boundary of domain in 1D
%   vdl : Dirichlet boundary condition of u
%   vnr : Newmann boundary condition of u
%   samples : # of ftn values used in drawing a graph
%   P: Order of approximation
%   viz: tag to turn on drawing
%
%   spdiffsineNB(-1, 1, 30, 5, 0, 0, 1);
%   spdiffsineNB(0, 2, 30, 5, 1, -1, 1);
%   spdiffsineNB(-50, 10, 80, 5, -40, 20, 1);
%   spdiffsineNB(-0.2, 0.4, 80, 5, 1, 2, 1);
%   spdiffsineNB(-.75, .75, 20, 5, 999, 0, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L, fvector, max_error] = spdiffsineNB(left, right, samples, P, vdl, vnr, viz)

% Standard parametrization: = [a, b]
% 
% d^2/dx^2 u = f where u(-1) = 0 = u(1) and f = sin(pi*x)

% Pseudo Jacobi Polynomial is defined as

%                          (1 - xi)/2    if p == 0
% pi[p](x) -> psi[p](xi) = (xi + 1)/2    if p == 1
%                          (1 - xi)/2 * (xi + 1)/2 *...
%                         P (1,1)[p-2](xi), if 2 <= p <= P
% , where xi is in [-1, 1].

% where x = XI(xi) = (a+b)/2 + (b-a)/2 * xi, where xi is in [-1, 1]
%  then, x = XI(xi) is in [a, b] where XI(-1) = a, XI(+1) = b.

% Since we need Integral[a, b](psi'[p](x)*(psi'[q](x)*dx, we change this integral on [a,b]
% to that on [-1, 1]. Then we can use Gaussian  quadrature for evaluation of the integral.

% Elemental Laplacian matrix is defined by L[p][q] = <psi_prime(p), psi_prime(q)>
% Since in PsdJcb1, we have almost diagonal except (0,P), (P,0)

% alpha, beta in Jacobi polynomials 
  alpha = 1;
  beta = 1;

% Define internal coefficients to be solved as
  u_hat = zeros(P-1, 1);

% Define Laplacian matrix L
  L = zeros(P-1, P-1);

% For RHS, define fvector[P-1] as zero
  fvector = zeros(P-1, 1);

% Define z,w for gaussian quadrature
  udegree = 2 * P; % p+1-1 + p+1-1 = 2p
  [dgq, dgrq, dglq] = QuadratureDegree(udegree);

% Gauss-Lobetto Quadrature->we take 'dglq' as a degree of corresponding Jacobi Polynomial
  qalpha = 0; qbeta = 0;
  [z, w] = JacobiGLZW(dglq, qalpha, qbeta);
  xiz = (right-left)/2 * z + (left+right)/2;

% Reparametrization variables
  Jn = (right-left)/2;
  
% Computing L and fvector
  for p = 2:P
    for q = 2:P
      
      mjcbvec1 = ModifiedJacobiPolyDerivative(p, z, alpha, beta);
      mjcbvec2 = ModifiedJacobiPolyDerivative(q, z, alpha, beta);
      % sum up for innerproduct w.uz - LHS
      L(p-1, q-1) = - w' * (mjcbvec1 .* mjcbvec2) / Jn;
    end

    % For RHS,
    fz = sin(pi*xiz) .* ModifiedJacobiPoly(p, z, alpha, beta);
    fvector(p-1, 1) = Jn * w' * fz;

  end

  uh_hat = inv(L) * fvector
    
  % Computing Neumann Boundary conditions
  fz = sin(pi*xiz) .* ModifiedJacobiPoly(1, z, alpha, beta);
  fvec = Jn * w' * fz;
    
  if vdl == 999
    bdl = -1/(pi^2) * sin(pi * left);
  else
    bdl = vdl;
  end
    
  u_hat0 = bdl;
  u_hat1 = u_hat0 + (left - right) * (fvec - vnr);
  u_hat = [u_hat0; u_hat1; uh_hat];

% For given force function f,    
% we have order of approximation P
%            h_hat vector[1,...,P-1]
%            psi[1,...,P-1](xi)
%            As an approximate solution of 
%                   arbitrary parametrization: = [left, right]
%                   d^2/dx^2 u = f 

% Plotting app.sol u on [left, right] and 
% compare with analytic solution -1/pi^2*sin(pi*xi)

max_error = 999;

if samples ~= 0
  N = samples; % Number of subintervals in [left, right]
   
  h = (right - left)/N;

  x = [left:h:right]';
  sx = 2/(right-left) * (x - 0.5*(left + right));
  
  MjcbPoly = ModifiedJacobiPoly(0, sx, alpha, beta);  % Starts from degree 0
  for l = 1:P   %Degree of Jacobi polynomial in sxi
    MjcbPolyDeg = ModifiedJacobiPoly(l, sx, alpha, beta);
    MjcbPoly = [MjcbPoly, MjcbPolyDeg];
  end
  
  u_xi = MjcbPoly * u_hat;

  for s = 1:P-1
    for t = 1:P-1
      elt = L(s,t);
      if abs(elt) < 1e-13
        L(s,t) = 0;
      end
    end
  end


  % Analytic solution considering boundary condition
  v = zeros(N+1, 1);

  tmp1 = 1/pi * cos(pi*right);
  tmp2 = 1/(pi^2) * sin(pi*left);
  
  C = vnr + tmp1;
  D = bdl - left * vnr + tmp2 - left * tmp1;

  v = -1/(pi^2) * sin(pi * x) + C * x + D;

  max_error = max(abs(u_xi - v))
    
  if viz ~= 0

  figure(1);
   subplot(3,1,1);
    plot(x, u_xi, '.-','markersize',13);
    grid on, title('Approximation by Spectral Method');
    %axis([left right -1.2 1.2]);
% 
%     
  subplot(3,1,2);
    plot(x, v, '.-','markersize',13);
    grid on, title('Analytic Solution');
    %axis([left right -0.2 0.2]);
    text(0, 0, num2str(max_error));
 
  subplot(3,1,3);
    spy(L);

  end
  
end

return  