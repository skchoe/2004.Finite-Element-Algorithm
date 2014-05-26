%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffstiffness.m Spectral Element Solver on 1D
%
%  example commands
%  left = 0; right = 200; range of whole domain
%  dxvec = [30, 30, 40, 40, 60]'; Vector showing length of element
%  pvec = [7, 3, 5, 9, 10]'; Vector showing orders on each element
%  vBL = 0, vBR = 0: Boundary conditions are zero
%  sigleft, sigright = Tensor values at each boundary
function [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(left, right, dxvec, pvec, vBL, vBR, kFun, fFun)
  
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

    % Initialize of solution domain, analytic solution image, numerical
    % solution Image
    xvec = [];
    yvec_apx = [];

    % Finding a table of quadrature points/weights for each element%%%%%%%%%
    % And Basis value table for each element%%%%%%%%%%%%%%%%%%%%
    qalpha = 0; qbeta = 0;  % These are (\alpha, \beta) for quadrature point calculation 
    Pmx = max(pvec(:,1));
    [dgq, dgrq, dglq] = QuadratureDegree(3 * Pmx);
    [z, w] = JacobiGLZW(dglq, qalpha, qbeta);
    Nzmx = size(z, 1);
    
    % Index Array for Number of Quadrature Points
    idxQpts = zeros(N, 1);
    
    % Quadrature Points /Weights
    zwTable = zeros(Nzmx, 2, N);
    [zwTable, idxQpts] = setup_zwTable(qalpha, qbeta, N, pvec, idxQpts, zwTable);
    
    % Basis Vaue Table Decleration
    basisTable = zeros(Nzmx, Pmx+1, N);
    basisTable = setup_basisTable(alpha, beta, N, idxQpts, pvec, zwTable, basisTable);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % for each element el = 1~N, compute global stiffness matrix Lg and RHS fg
    Pn = 0;  % Order of each element -> Order of the last element
    for el = 1:N

        Pn = pvec(el, 1); % order for each element el

        % Making eleft/eright
        eleft = left;
        eright = right;
        for i = 1:el-1
            eleft = eleft + dxvec(i);
        end
        eright = eleft + dxvec(el);
        if eright > right
            eright = right;
        end

        Nqpts = idxQpts(el, 1);
        z = zwTable(1:Nqpts, 1, el);
        w = zwTable(1:Nqpts, 2, el);
        xiz = (eright-eleft)/2 * z + (eleft+eright)/2;

        % Global x vector definition
        xvec = [xvec; xiz];
    
        % Jacobian in Reparametrization variables
        Jn = abs(eright-eleft)/2;

        [L, fvector] = spLocalMats1D(Pn, alpha, beta, z, w, kFun, fFun, xiz, Jn);

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
    [Lg, Fg] = ApplyBdy2Matrix(left, right, Ndof, Pn, Lg, Fg, vBL, vBR, kFun);
  
    L = inv(Lg);
    % global coefficient ug_hat
    ug_hat = inv(Lg) * Fg;

    for iel = 1:N
        p_i = pvec(iel);
        coeloc = coeffglobal2local(ug_hat, map, pvec, iel);
        
        Nqpts = idxQpts(iel, 1);
        bt = basisTable(1:Nqpts,1:p_i+1,iel);
        uapx = coeloc*bt';
        yvec_apx = [yvec_apx;uapx'];

    end

    % final plotting evaluating the linear expansion of solution with respect to basis functions
    figno = 33;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(figno);

        % Numerical solution considering boundary condition
        plot(xvec, yvec_apx, '.-','markersize',5);
        grid on, title('Numerical Solution');


return

function [Lg, Fg] = ApplyBdy2Matrix(left, right, Ndof, Pn, Lg, Fg, vBL, vBR, kFun)

    idEnd = Ndof;   
    kR = feval(kFun, right)
    kL = feval(kFun, left)

    vecl = Lg(:, 1);        vecl(1, 1) = -1;
    vecr = zeros(Ndof, 1);  vecr(idEnd, 1) = 1;

    % Insert special value on Lg and Fg for boundary condition
    Lg(1, :) = 0; Lg(:, 1) = 0;   Lg(1, 1) = 1; 
  
    Fg(1, 1) = 0;

    RT = kR * vecr * vBR;
    LT = vecl * vBL;
  
    Fg = Fg + RT - LT;

return


function [zwTable, idxQpts] = setup_zwTable(qalpha, qbeta, N, pvec, idxQpts, zwTable)

    for ei = 1:N
        p_i = pvec(ei, 1);
        
        % Define z, w for gaussian quadrature
        udegree = 3 * p_i; % degree of Kappa + p+1-1 + p+1-1 = 3p approximately
        [dgq, dgrq, dglq] = QuadratureDegree(udegree);
        
        % Gauss-Lobetto Quadrature->we take 'dglq' as a degree of corresponding Jacobi Polynomial
        
        [z, w] = JacobiGLZW(dglq, qalpha, qbeta);
        sz = size(z, 1);
        
        idxQpts(ei, 1) = sz;
        
        for j = 1:sz
            zwTable(j, 1, ei) = z(j, 1);
            zwTable(j, 2, ei) = w(j, 1);
        end
    end
return

function basisTable = setup_basisTable(alpha, beta, N, idxQpts, pvec, zwTable, basisTable)
    for ei = 1:N
        p_i = pvec(ei, 1);
        z_i = idxQpts(ei, 1);
        
        for ip = 0:p_i
            for iz = 1:z_i
                z = zwTable(iz, 1, ei);
                basisTable(iz, ip+1, ei) = ModifiedJacobiPoly(ip, p_i, z, alpha, beta);
            end
        end
    end

return

function coeloc = coeffglobal2local(ug_hat, map, pvec, iel)
    odr = pvec(iel, 1);
    coeloc = zeros(1, odr+1);
    
    idxvec = map(iel, :);
    
    for i = 1 : odr+1
        k = idxvec(1, i);
        coeloc(1, i) = ug_hat(k, 1);
    end
return