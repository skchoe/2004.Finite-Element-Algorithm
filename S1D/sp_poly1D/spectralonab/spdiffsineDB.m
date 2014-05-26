%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffsine.m Spectral Element Solver on 1D
%
%   max_error = spdiffsineDB(left, right, samples, P, vdl, vdr, viz)
%   left, right : boundary of domain in 1D
%   vdl, vdr : Dirichlet boundary conditions of u
%   samples : # of ftn values used in drawing a graph
%   P: Order of approximation
%   viz: tag to turn on drawing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Standard parametrization: = [a, b]
% 
% d^2/dx^2 u = f where u(-1) = 0 = u(1) and f = sin(pi*x)
%
% Pseudo Jacobi Polynomial is defined as
%
%                          (1 - xi)/2    if p == 0
% pi[p](x) -> psi[p](xi) = (1 + xi)/2    if p == 1
%                          (1 - xi)/2 * (xi + 1)/2 *...
%                         P (1,1)[p-2](xi), if 2 <= p <= P
% , where xi is in [-1, 1].
%
% where x = XI(xi) = (a+b)/2 + (b-a)/2 * xi, where xi is in [-1, 1]
%  then, x = XI(xi) is in [a, b] where XI(-1) = a, XI(+1) = b.
%
% Since we need Integral[a, b](psi'[p](x)*(psi'[q](x)*dx, we change this integral on [a,b]
% to that on [-1, 1]. Then we can use Gaussian  quadrature for evaluation of the integral.
%
% Elemental Laplacian matrix is defined by L[p][q] = <psi_prime(p), psi_prime(q)>
% Since in PsdJcb1, we have almost diagonal except (0,P), (P,0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% left = -1; right = 1; samples = 90; P = 3; vdl = 0; vdr = 0; viz = 1;
% left = 0; right = 1; samples = 30; P = 8; vdl = 0; vdr = 0; viz = 1;
% left = -5; right = 1; samples = 80; P = 5; vdl = -40; vdr = 20; viz = 1;
%left = -.5; right = 0.0; samples = 80; P = 5; vdl = 0.1013; vdr = - 0.1013 ; viz = 1;
%left = -.75; right = .75; samples = 20; P = 3; vdl = 0.0716; vdr = - 0.0716 ; viz = 1;
 left = -0.2; right = 0.4; samples = 80; P = 5; vdl = 1; vdr = 2; viz = 1;
% spdiffsineDB(-1, 1, 90, 3, 0, 0, 1);
%function [L, fvector, error] = spdiffsineDB(left, right, samples, P, vdl, vdr, viz)


% vdl = -1/(pi^2) * sin(pi*left);
% vdr = -1/(pi^2) * sin(pi*right);
% 
% % alpha, beta in Jacobi polynomials 
%   alpha = 1;
%   beta = 1;
% 
% % Define internal coefficients to be solved as
%   u_hat = zeros(P+1, 1);
% 
% % Define Laplacian matrix L
%   L = zeros(P+1, P+1);
% 
% % For RHS, define fvector[P-1] as zero
%   fvector = zeros(P+1, 1);
% 
% % Define z,w for gaussian quadrature
%   udegree = 2 * P; % p+1-1 + p+1-1 = 2p
%   [dgq, dgrq, dglq] = QuadratureDegree(udegree);
% 
% % Gauss-Lobetto Quadrature->we take 'dglq' as a degree of corresponding Jacobi Polynomial
%   qalpha = 0; qbeta = 0;
%   [z, w] = JacobiGLZW(dglq, qalpha, qbeta);
%   xiz = (right-left)/2 * z + (left+right)/2;
% 
%   fzs = sin(pi*xiz);
% 
%   % Jacobian in Reparametrization variables
%   Jn = (right-left)/2;
% 
%   [L, fvector] = spLocalMats1D(P, alpha, beta, z, w, fzs, Jn);
%   
%   vl = L(:, 1); vl(1, 1) = -1;   vl(2, 1) = 0;
%   vr = L(:, 2); vr(1, 1) = 0;   vr(2, 1) = -1;
%   
%   L(1, :) = 0;  L(:, 1) = 0;    L(1, 1) = 1;
%   L(2, :) = 0;  L(:, 2) = 0;    L(2, 2) = 1;
%   
%   fvector(1, 1) = 0;    fvector(2, 1) = 0;
%   
% 
%   F = fvector - vl*vdl - vr*vdr;
% 
%   u_hat = inv(L) * F;
% 
%   
% Plotting app.sol u on [left, right] and
% compare with analytic solution -1/pi^2*sin(pi*xi)
% 
% max_error = 999;
% if samples ~= 0
%   N = samples; % Number of subintervals in [left, right]
% 
%   [xvec, yvec] = spEvalSpan1D(P, alpha, beta, left, right, N, u_hat);
%   
%   L = spSetSpy(L, 10e-14, 0, 0);
% 
%   % Analytic solution considering boundary condition
%   yavec = zeros(N+1, 1);
%   tmpL = -1/(pi^2) * sin(pi*left);
%   tmpR = -1/(pi^2) * sin(pi*right);
%   C = (vdr - vdl + tmpL - tmpR)/(right - left);
%   D = vdr - tmpR - C*right;
%   
%   yavec = -1/(pi^2) * sin(pi * xvec) + C * xvec + D;
% 
%   max_error = max(abs(yvec - yavec))
%     
%   if viz ~= 0
% 
%   figure(1);
%   subplot(3, 1, 1);
%     plot(xvec, yvec, '.-','markersize',13);
%     grid on, title('Approximation by Spectral Method');
%     %axis([left right -1.2 1.2]);
% 
%     
%   subplot(3, 1, 2);
%     plot(xvec, yavec, '.-','markersize',13);
%     grid on, title('Analytic Solution');
%     %axis([left right -0.2 0.2]);
%     text(0, 0, num2str(max_error));
%  
%   subplot(3, 1, 3);
%     spy(L);
% 
%   end
%   
% end
return