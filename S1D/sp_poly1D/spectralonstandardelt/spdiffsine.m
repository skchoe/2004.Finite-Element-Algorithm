%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffsine.m Spectral Element Solver on 1D
%
%
% spdiffsine(P, viz, vl, vr, samples)
% P: Order of polynomials
% vl, vr: boundary values of left and right repectively.
% by default, u(-1)=1, u(1) =-1.
% example
% spdiffsine(5, 1, 1, -1, 50);
% spdiffsine(5, 1, -3, 5, 50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function max_error = spdiffsine(P, viz, vl, vr, samples)

% Standard parametrization: = [-1, 1]
right = 1;
left = -1;

% d^2/dx^2 u = f where u(-2) = 0 = u(3) and f = sin(pi*x)

% Pseudo Jacobi1 is defined as

%                           (1 - xi)/2    if p == 0
% pi[p](xi) -> psi[p](xi) = (1 + xi)/2    if p == 1
%                           (1 - xi)/2 * (1 + xi)/2 * P(1,1)[p-2](xi), if 2 <= p <= P

% Elemental Laplacian matrix is defined by L[p][q] = <psi_prime(p), psi_prime(q)>
% Since in PsdJcb1, we have almost diagonal except (0,P), (P,0)

% We use gauss-Lobatto quadrature for L[p][q] where 1 < p, q <= P
% By boundary condition, u_hat(0) = u_hat(1) = 0

% Define the maximum order P of the element e
  alpha = 1;
  beta = 1;

% Define internal coefficients to be solved as
  uh_hat = zeros(P-1, 1);

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
  
  fz1 = sin(pi*z) .* ModifiedJacobiPoly(0, z, alpha, beta);
  fvec0 = w'*fz1
  fz2 = sin(pi*z) .* ModifiedJacobiPoly(1, z, alpha, beta);
  fvec1 = w'*fz2
  
  mjcbvec0 = ModifiedJacobiPolyDerivative(0, z, alpha, beta);
  mjcbvec1 = ModifiedJacobiPolyDerivative(1, z, alpha, beta);
  L00 = -w' * (mjcbvec0 .* mjcbvec0)
  L10 = -w' * (mjcbvec1 .* mjcbvec0)
  L01 = -w' * (mjcbvec0 .* mjcbvec1)
  L11 = -w' * (mjcbvec1 .* mjcbvec1)
  
  
  % Note that the matrix L is invertible
    uh_hat = inv(L) * fvector;

    u_hat = [vl, vr, uh_hat']';
    
% For given force function f,    
% we have order of approximation P
%            h_hat vector[1,...,P-1]
%            psi[1,...,P-1](xi)
%            As an approximate solution of 
%                   Standard parametrization: = [-1, 1]
%                   d^2/dx^2 u = f where u(-1) = 0 = u(1)

% Plotting app.sol u on [-1, 1] and 
% compare with analytic solution -1/pi^2*sin(pi*xi)

  N = samples; % Number of subintervals in [-1, 1]
    
  h = (right - left)/N;
  xi = [left:h:right]'; % Samples in x-axis

  MjcbPoly = ModifiedJacobiPoly(0, xi, alpha, beta);  % Starts from degree 0
  for l = 1:P   %Degree of Jacobi polynomial in sxi
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
 
  % case which oundaries are 0, 0.
  tag = 1;
  % case which oundaries are 1, -1.
  if vl==1 & vr == -1
    tag = 0;
  elseif vl==-3 & vr == 5
    tag = 2;
  end
  
  v = zeros(N+1, 1);
  for n = 1:N+1
    if tag == 1
        v(n, 1) = -1/(pi^2) * sin(pi * xi(n,1));
    elseif tag == 0
        v(n, 1) = -1/(pi^2) * sin(pi * xi(n,1)) - xi(n,1);
    elseif tag == 2
        v(n, 1) = -1/(pi^2) * sin(pi * xi(n,1)) + 4 * xi(n,1) + 1;
    end    
  end

  max_error = max(abs(u_xi - v));
    
  if viz ~= 0

  figure(1);
  subplot(2,1,1);
    plot(xi, u_xi, '.-','markersize',13);
    grid on, title('Approximation by Spectral Method');
    %axis([-1, 1, -.2, .2]);

    
  subplot(2,1,2);
    plot(xi, v, '.-','markersize',13);
    grid on, title('Analytic Solution');
    %axis([-1, 1, -.2, .2]);

    text(0, 0, num2str(max_error));
%  
%   subplot(2,1,3);
% 
%   spy(L);

  end
return  