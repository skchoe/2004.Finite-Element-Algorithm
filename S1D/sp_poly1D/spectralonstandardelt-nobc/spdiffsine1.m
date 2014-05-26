%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffsine.m Spectral Element Solver on 1D
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function max_error = spdiffsine1(P, viz)

% Standard parametrization: = [-1, 1]
left = -1;
right = 1;

% d^2/dx^2 u = f where u(-2) = 0 = u(3) and f = sin(pi*x)

% Pseudo Jacobi1 is defined as

%                           (1 - xi)/2    if p == 0
% pi[p](xi) -> psi[p](xi) = (1 + xi)/2    if p == 1
%                           (1 - xi)/2 * (1 + xi)/2 * P(1,1)[p-1](xi), if 2 <= p <= P

% Elemental Laplacian matrix is defined by L[p][q] = <psi_prime(p), psi_prime(q)>
% Since in PsdJcb1, we have almost diagonal except (0,P), (P,0)

% We use gauss-Lobatto quadrature for L[p][q] where 1 < p, q <= P
% By boundary condition, u_hat(0) = u_hat(1) = 0

% Define the maximum order P of the element e
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
  
% Computing L and fvector
  for p = 2:P
    for q=2:P,
      
      mjcbvec1 = ModifiedJacobiPolyDerivative(p, z, alpha, beta);
      mjcbvec2 = ModifiedJacobiPolyDerivative(q, z, alpha, beta);
      % sum up for innerproduct w.uz - LHS
      L(p-1, q-1) = - w' * (mjcbvec1 .* mjcbvec2);
    end

    % For RHS,
    fz = sin(pi*z) .* ModifiedJacobiPoly(p, z, alpha, beta);
    fvector(p-1, 1) = w' * fz;

  end
  L

  % Note that the matrix L is invertible
    u_hat = inv(L) * fvector;

% % % % %   % Testing point by point
% % % % %   x = 0.3;
% % % % %   vx = -1 / pi^2 * sin(pi*x);
% % % % %   ux = 0;
% % % % %   for i = 2:P
% % % % %     ux = ux + u_hat(i-1, 1)*ModifiedJacobiPoly(i, x, alpha, beta);
% % % % %   end
% % % % %   Error = abs(ux-vx)


% For given force function f,    
% we have order of approximation P
%            h_hat vector[1,...,P-1]
%            psi[1,...,P-1](xi)
%            As an approximate solution of 
%                   Standard parametrization: = [-1, 1]
%                   d^2/dx^2 u = f where u(-1) = 0 = u(1)

% Plotting app.sol u on [-1, 1] and 
% compare with analytic solution -1/pi^2*sin(pi*xi)

  N = 30; % Number of subintervals in [-1, 1]
  xi = zeros(N+1, 1);
%   u_xi = zeros(N+1, 1);
%   
%   % Boundary condition u(right) = 0 = u(left)
%   u_xi(1, 1) = 0;
%   u_xi(N+1, 1) = 0;
    
  h = (right - left)/N;
  
  xi = [left:h:right]';

  MjcbPoly = ModifiedJacobiPoly(2, xi, alpha, beta);  % Starts from degree 0
  for l = 3:P   %Degree of Jacobi polynomial in sxi
    MjcbPolyDeg = ModifiedJacobiPoly(l, xi, alpha, beta);
    MjcbPoly = [MjcbPoly, MjcbPolyDeg];
  end
  
  u_xi = MjcbPoly * u_hat;

  % compute L for spy();
  for s = 1:P-1
    for t = 1:P-1
      elt = L(s,t);
      if abs(elt) < 1e-13
        L(s,t) = 0;
      end
    end
  end
    
  v = zeros(N+1, 1);
  for n = 1:N+1
    v(n, 1) = -1/(pi^2) * sin(pi * xi(n,1));
  end

  max_error = max(abs(u_xi - v));
    
  if viz ~= 0

  figure(1);
  subplot(3,1,1);
    plot(xi, u_xi, '.-','markersize',13);
    grid on, title('Approximation by Spectral Method');
    %axis([left, right, -1.2, 1.2]);

    
  subplot(3,1,2);
    plot(xi, v, '.-','markersize',13);
    grid on, title('Analytic Solution');
    %axis([left, right, -1.2, 1.2]);

    text(0, 0, num2str(max_error));
 
  subplot(3,1,3);

  spy(L);

  end
return  