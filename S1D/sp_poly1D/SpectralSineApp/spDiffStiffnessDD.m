%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffstiffness.m Spectral Element Solver on 1D
%
% example commands
%  left = 0; right = 200; range of whole domain
%  dxvec = [30, 30, 40, 40, 60]'; Vector showing length of element
%  pvec = [7, 3, 5, 9, 10]'; Vector showing orders on each element
%  vdl = 0, vdr = 0: Boundary conditions are zero
%
 function max_error = spDiffStiffnessDD(left, right, dxvec, pvec, vdl, vdr, numsampleintervals)
%  left = -.75; right = .75; dxvec = [1.5]'; pvec = [5]'; vdl = 999; vdr = 999; numsampleintervals = 20;
%  left = -5; right = 1; dxvec = [6]'; pvec = [1]'; vdl = 999; vdr = 999; numsampleintervals = 40;
%  left = -.5; right = 1.0; dxvec = [.5, 1.0]'; pvec = [10,15]'; vdl = 999; vdr = 999; numsampleintervals = 100; 
%    left = -.5; right = 1.0; dxvec = [.25, .25, 1.0]'; pvec = [1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 10; 
%    left = -.5; right = 1.0; dxvec = [.25, .25, .5, .5]'; pvec = [1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 10; 
%    left = -.5; right = 1.0; dxvec = [.25, .15, .1, .1, .4, .5]'; pvec = [1, 1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 10; 
%    left = -.5; right = 1.0; dxvec = [.25, .15, .1, .1, .3, .3, .3]'; pvec = [1, 1, 1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 10; 
%    left = -.5; right = 1.0; dxvec = [.1, .1, .1, .1, .1, .3, .3, .3]'; pvec = [1, 1, 1, 1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 10; 
%    left = -.5; right = 1.0; dxvec = [.1, .1, .1, .1, .1, .1, .2, .1, .05, .05, .1, .3]'; pvec = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 10; 
%  left = 0; right = 2; dxvec = [.5,.5,.5, .5]';pvec = [1, 1,1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 10; 
%  left = 0; right = 2; dxvec = [.3, .3, .4, .4, .6]'; pvec = [7, 3, 5, 9, 10]'; vdl = 999; vdr = 999; numsampleintervals = 100; 
%  left = 0; right = 10; dxvec = [1, 2, 3, 1, 3]'; pvec = [7, 3, 5, 9, 10]'; vdl = 999; vdr = 999; numsampleintervals = 100; 
%  left = 0; right = 10; dxvec = [1.5, 1.5, 2.8, 1.2, 3]'; pvec = [1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 100; 
%  left = 0; right = 20; dxvec = [3, 3, 4, 4, 6]'; pvec = [7, 3, 5, 9, 10]'; vdl = 999; vdr = 999; numsampleintervals = 100; 
%  left = 0; right = 20; dxvec = [3, 3, 4, 4, 6]'; pvec = [1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 100; 
%   left = -2.5; right = 2.5; dxvec = [1, 1, 1, 1, 1]'; pvec = [9, 9, 9, 9, 9]'; vdl = 999; vdr = 999; numsampleintervals = 100;
%   left = -2.5; right = 2.5; dxvec = [1, 1, 1, 1, 1]'; pvec = [8, 8, 8, 8, 8]'; vdl = 999; vdr = 999; numsampleintervals = 20;
%   left = -2.5; right = 2.5; dxvec = [1, 1, 1, 1, 1]'; pvec = [5, 5, 5, 5, 5]'; vdl = 999; vdr = 999; numsampleintervals = 100;
%   left = -2.5; right = 2.5; dxvec = [1, 1, 1, 1, 1]'; pvec = [4, 4, 4, 4, 4]'; vdl = 999; vdr = 999; numsampleintervals = 100;
%   left = -2.5; right = 2.5; dxvec = [1, 1, 1, 1, 1]'; pvec = [3, 3, 3, 3, 3]'; vdl = 999; vdr = 999; numsampleintervals = 100;
%   left = -2.5; right = 2.5; dxvec = [1, 1, 1, 1, 1]'; pvec = [2, 2, 2, 2, 2]'; vdl = 999; vdr = 999; numsampleintervals = 100;
%   left = -2.5; right = 2.5; dxvec = [1, 1, 1, 1, 1]'; pvec = [1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 100;
%   left = -2.5; right = 2.5; dxvec = [.5, .5, .5, .5, .5, .5, .5, .5, .5, .5]'; pvec = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]'; vdl = 999; vdr = 999; numsampleintervals = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   
% alpha, beta in Jacobi polynomials
  alpha = 1;
  beta = 1;

  % error handling
  sdxvec = size(dxvec, 1);
  spvec = size(pvec, 1);

  % number of elements
  N = size(dxvec, 1);

  % define map from local to global stiffness marix using modified Jacobi polynomials with
  % linear mode on indices 1, 2. From 3, bubble modes are assigned.
  map = spCalcMap(N, pvec);

  % determine maximum dimension of stiffness matrix Lg
  Ndof = pvec(1, 1)+1;
  for i = 2:N
    Ndof = Ndof + pvec(i, 1);
  end

  % define stiffness matrix Lg and RHS Fg
  Lg = zeros(Ndof, Ndof);
  Fg = zeros(Ndof, 1);

  % for each element el = 1~N, compute global stiffness matrix Lg and RHS fg
  Pn = 0;  % Order of each element -> Order of the last element
  for el = 1:N

    Pn = pvec(el, 1); % order for each element el

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eleft = left;
    eright = right;

    for i = 1:el-1
      eleft = eleft + dxvec(i);
    end

    eright = eleft + dxvec(el);

    if eright > right
      eright = right;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define z, w for gaussian quadrature
    udegree = 2 * Pn; % p+1-1 + p+1-1 = 2p
    [dgq, dgrq, dglq] = QuadratureDegree(udegree);

    % Gauss-Lobetto Quadrature->we take 'dglq' as a degree of corresponding Jacobi Polynomial
    qalpha = 0; qbeta = 0;
    [z, w] = JacobiGLZW(dglq, qalpha, qbeta);
    xiz = (eright-eleft)/2 * z + (eleft+eright)/2;

    fzs = sin(pi*xiz);

    % Jacobian in Reparametrization variables
    Jn = (eright-eleft)/2;

    [L, fvector] = spLocalMats1D(Pn, alpha, beta, z, w, fzs, Jn);

    % composing global stiffness matrix Lg(pg, qg)
    for p = 1:Pn + 1
        
      pg = map(el, p);
      % RHS term processing
      Fg(pg, 1) = Fg(pg, 1) + fvector(p, 1);

      % Lg(pg, qg) processing
      for q = 1:Pn + 1
        qg = map(el, q);
        Lg(pg, qg) = Lg(pg, qg) + L(p, q);
      end
    end

  end
  

  % Find index where right-most mode reside

  if Ndof == 2
    idEnd = 2; 
  else
    idEnd = Ndof - Pn + 1; 
  end

  vecl = Lg(:, 1);  vecl(1, 1) = -1;    vecl(idEnd, 1) = 0;

  vecr = Lg(:, idEnd);  vecr(1, 1) = 0; vecr(idEnd, 1) = -1;

  % Insert special value on Lg and Fg for boundary condition
  Lg(1, :) = 0; Lg(:, 1) = 0;   Lg(1, 1) = 1; 
  Lg(idEnd, :) = 0; Lg(:, idEnd) = 0; Lg(idEnd, idEnd) = 1;

  Fg(1, 1) = 0; Fg(idEnd, 1) = 0;

  if vdl == 999
    bdleft = -1/(pi^2) * sin(pi * left);
  else
    bdleft = vdl;
  end
  
  if vdr == 999
    bdright = -1/(pi^2) * sin(pi * right);
  else 
    beright = vdr;
  end
  
  
  A = vecl * bdleft;
  B = vecr * bdright;

  FG = Fg - A - B;

  % global coefficient ug_hat
  ug_hat = inv(Lg) * FG;

  % final plotting evaluating the linear expansion of solution with respect to basis functions
  figno = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  max_error = spPlotResultDD(numsampleintervals, left, right, alpha, beta, dxvec, pvec, Lg, ug_hat, map, figno);


