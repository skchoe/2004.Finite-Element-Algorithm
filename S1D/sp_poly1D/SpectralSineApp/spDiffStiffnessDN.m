%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffstiffness.m Spectral Element Solver on 1D
%
% example commands
%  left = 0; right = 200; range of whole domain
%  dxvec = [30, 30, 40, 40, 60]'; Vector showing length of element
%  pvec = [7, 3, 5, 9, 10]'; Vector showing orders on each element
%  vdl = 0, vnr = 0: Boundary conditions are zero
%
function max_error = spDiffStiffnessDN(left, right, dxvec, pvec, vdl, vnr, numsampleintervals)
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

    [uvec, fzs, vLDB, vRNB] = spEvalUserFunctions(left, right, xiz, size(xiz,1), vdl, vnr);
    
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

  if vdl == 999
    bdleft = vLDB;
  else
    bdleft = vdl;
  end
  if vnr == 999
    bdright = vRNB;
  else
    bdright = vnr;
  end

  LT = vecl * bdleft;
  RT = vecr * bdright;

  FG = Fg - RT - LT;

  % global coefficient ug_hat
  ug_hat = inv(Lg) * FG;

  % final plotting evaluating the linear expansion of solution with respect to basis functions
  figno = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  max_error = spPlotResultDN(numsampleintervals, left, right, alpha, beta, dxvec, pvec, bdleft, bdright, Lg, ug_hat, map, figno);
  
  