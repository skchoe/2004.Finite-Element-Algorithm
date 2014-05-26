%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spStiffCurveDN.m Spectral Element Solver on 1D
%
% example commands
%  left = 0; right = 200; range of whole domain
%  dxvec = [30, 30, 40, 40, 60]'; Vector showing length of element
%  pvec = [7, 3, 5, 9, 10]'; Vector showing orders on each element
%  vdl = 0, vnr = 0: Boundary conditions are zero
%
function max_error = spStiffCurveDN(left, right, prob_order, dxvec, pvec, vdl, vnr, numsampleintervals)
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [1, 1, 12, 1, 1]'; vdl = 0; vnr = 0; numsampleintervals = 300; %ERR(order = 5) = 0.0255
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [2, 2, 2, 2, 2]'; vdl = 0; vnr = 0; numsampleintervals = 300;  %ERR(order = 5) = 0.0019
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [3, 3, 3, 3, 3]'; vdl = 0; vnr = 0; numsampleintervals = 300;  %ERR(order = 5) = 2.4000e-004
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [4, 4, 4, 4, 4]'; vdl = 0; vnr = 0; numsampleintervals = 300;  %ERR(order = 5) = 5.6437e-006
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [5, 5, 5, 5, 5]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 5) = 5.7732e-015
%prob_order = 5;
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [5, 5, 5, 5, 5]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 2.6667e-006
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [7, 7, 7, 7, 7]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.2412e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [8, 8, 8, 8, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.0036e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [9, 9, 9, 9, 9]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 9.9920e-014
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [10, 10, 10, 10, 10]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.0636e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [11, 11, 11, 11, 11]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.0025e-013
   
   
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [5, 5, 5, 5, 5]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15) = 4.3229e-006
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [7, 7, 7, 7, 7]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15) = 2.2388e-007
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [8, 8, 8, 8, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15) = 8.6603e-008
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [9, 9, 9, 9, 9]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15) = 6.7875e-008
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [11, 11, 11, 11, 11]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15) = 6.6704e-008
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [13, 13, 13, 13, 13]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15)=  6.6711e-008
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [15, 15, 15, 15, 15]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15)=  6.6767e-008
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [8, 21, 21, 21, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15)=  takes everlasting
% 
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [8, 8, 15, 15, 15, 15, 15, 15, 8, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 19) =  6.1087e-004
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [8, 8, 15, 15, 15, 15, 15, 15, 8, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 15) =  6.6703e-008
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [8, 8, 15, 15, 15, 15, 15, 15, 8, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 10) =  2.4680e-011
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [9, 9, 16, 16, 16, 16, 16, 16, 9, 9]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 10) =  1.7741e-011
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [10, 10, 19, 19, 19, 19, 19, 19, 10, 10]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 10) =  1.8108e-011
% 
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [8, 8, 15, 15, 15, 15, 15, 15, 8, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 9) = 6.4685e-012
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [9, 9, 16, 16, 16, 16, 16, 16, 9, 9]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 9) = 5.7636e-012

% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [5, 5, 6, 6, 6, 6, 6, 6, 5, 5]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 5) =  7.5495e-015
% left = 0; right = 1; dxvec = [.1, .1, .1, .1, .1,.1, .1, .1, .1, .1]'; pvec = [5, 5, 6, 6, 6, 6, 6, 6, 5, 5]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 5) =  

% Critical example for spectral method
%prob_order = 5;


[coefvec, exactsolvec] = spPoissCurve(prob_order, numsampleintervals);

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

    % RHS evaluation of which coefficients were computed on top.
    fzs = spEvalPolynomial(coefvec, prob_order-1, xiz, size(z, 1));

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

  vecl = Lg(:, 1);  vecl(1, 1) = -1;

  vecr = zeros(Ndof, 1);    vecr(idEnd, 1) = 1;

  % Insert special value on Lg and Fg for boundary condition
  Lg(1, :) = 0; Lg(:, 1) = 0;   Lg(1, 1) = 1; 
  
  % Set Fg's DBC part to zero, no changes for NBC on right side
  Fg(1, 1) = 0;

  % Setting BC's: Left==Dirichlet, Right==Neumann BC
  if vdl == 999  %% DBC is not specified
    bdleft = -1/(pi^2) * sin(pi * left);
  else
    bdleft = vdl;
  end
  
  bdright = vnr;

  LT = vecl * bdleft;
  RT = vecr * bdright;

  FG = Fg - RT - LT;

  % global coefficient ug_hat
  ug_hat = inv(Lg) * FG;

  % final plotting evaluating the linear expansion of solution with respect to basis functions
  figno = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  max_error = spPoissCurveDN(numsampleintervals, left, right, alpha, beta, dxvec, pvec, bdleft, bdright, Lg, ug_hat, map, figno, exactsolvec);