%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANNdiffstiffness.m Spectral Element Solver on 1D
%
% example commands
%  BdyType: DD or DN type
%  left = 0; right = 200; range of whole domain
%  dxvec = [30, 30, 40, 40, 60]'; Vector showing length of element
%  pvec = [7, 3, 5, 9, 10]'; Vector showing orders on each element
%  vBL = 0, vBR = 0: Boundary conditions are zero
%  Nth = 8 Number of divisions on \theta direction
%  Nr = Number of divisions on r dirction
function [xvec, yvec, zMat_apx, vBL, vBR] = ANN2DSolvePoisson(BdyType, left, right, dxvec, pvec, Nth, vecBL, vecBR, sigleft, sigright)

    % h: interval in \theta direction
    h = 2 * pi / Nth;
    yvec = [0:h:2*pi-h]';

    % On R direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = 1;  % alpha, beta in Jacobi polynomials
    beta = 1;
    Nr = size(dxvec, 1);        % number of elements
    
    % define map from local to global stiffness marix using modified Jacobi polynomials with
    map = ANNCalcMap(Nr, pvec);  % linear mode on indices 1, N. bubble modes 2 ~ N-1.

    Ndof = 1;        % determine maximum dimension of stiffness matrix Lg: Global Deg. Freedom
    for i = 1:Nr        Ndof = Ndof + pvec(i, 1);       end

    % Finding a table of quadrature points/weights for each element%%%%%%%%
    qalpha = 0; qbeta = 0;  % These are (\alpha, \beta) for quadrature point calculation 
    Pmx = max(pvec(:,1));
    [dgq, dgrq, dglq] = QuadratureDegree(2 * Pmx);
    [z, w] = JacobiGLZW(dglq, qalpha, qbeta);
    Nzmx = size(z, 1);  % Maximum length of z,w values among elements

    % Index Array for Number of Quadrature Points for each element
    % nqpoints: Number of quadrature points on each element
    nqpoints = zeros(Nr, 1);

    % Quadrature Points /Weights
    zwTable = zeros(Nzmx, 2, Nr);
    [zwTable, nqpoints] = setup_zwTable(qalpha, qbeta, Nr, pvec, nqpoints, zwTable);

    % Local Matrix Table - numMat(3) Matrices for each element
    mx = max(pvec) + 1;   % Maximum local order + 1 = Max deg. freedom of elements 
    numMat = 3;           % Indicate M1, M2, M3 matrices/
    TL = zeros(mx, mx, numMat, Nr); % Storage of M1, M2, M3 for each element/
    Fel = zeros(Nth, mx, Nr);       %Storage for \int r^2 f(r,\theta) \psi \exp dr dtheta

    [TL, Fel, xvec] = setup_LocMatTable(BdyType, left, right, alpha, beta, Nr, Nth, ...
                                  dxvec, pvec, zwTable, nqpoints, numMat, TL, Fel, vecBL, vecBR);

    %-----------------------------------------------------------------------------
    % Starting local matrix/rhs and global assembly for each \theta-direction element

    % Global Stiffness Matrix considering r, \theta
    Mg = zeros(Ndof*Nth, Ndof*Nth);
    
    % LHS global assembly   
    for k = 1 : Nth
        % Global Stiffness Matrix considering only r
        Mgloc = zeros(Ndof, Ndof);

        p_ = 0;
        for ei = 1:Nr
            p_ = pvec(ei, 1);
            M1 = TL(1:p_+1, 1:p_+1, 1, ei);
            M2 = TL(1:p_+1, 1:p_+1, 2, ei);
            M3 = TL(1:p_+1, 1:p_+1, 3, ei);

            % Defining appropriate K
            Hf = Nth/2;
            j = k-1;
            if j > Hf
                K = -Hf + (j-Hf);
            else
                K = j;
            end

            % Local matrix Mel
            Mel = M3 + M2 + K^2 * M1;

            for p = 1:p_ + 1
                pg = map(ei, p);
   
                for q = 1:p_ + 1
                    qg = map(ei, q);
                    Mgloc(pg, qg) = Mgloc(pg, qg) + Mel(p, q);
                end
    
            end

        end
        
        loc_beg = Ndof * (k - 1) + 1;
        loc_end = loc_beg + Ndof - 1;
        Mg(loc_beg:loc_end, loc_beg:loc_end) = Mgloc;
        
    end
    
    % RHS global assembly
    Fg = zeros(Ndof*Nth, 1);
    Fgloc = zeros(Ndof, Nth);
    for k = 1 : Nth
        
        p_ = 0;
        for ei = 1:Nr
            p_ = pvec(ei, 1);

            for p = 1:p_ + 1
                pg = map(ei, p);
   
                Fgloc(pg, k) = Fgloc(pg, k) + Fel(k, p, ei);
            end

        end
        
        loc_beg = Ndof * (k - 1) + 1;
        loc_end = loc_beg + Ndof - 1;
        Fg(loc_beg:loc_end, 1) = Fgloc(:, k);      
        
    end
    

    % IFFT/FFT on Boundary terms
    % Need FFT to get u_hat_{0,k}, To apply Neumann cond. need IFFT(vBR).

    vecBL = fft(vecBL / Nth);
    vecBR = fft(vecBR) / Nth;


%   Boundary Condition Processing
    [Mg, Fg] = applyBdy2GMat(BdyType, left, right, Ndof, pvec, Nth, Mg, Fg, vecBL, vecBR, sigleft, sigright);


    % global coefficient ug_hat
    ug_hat = inv(Mg) * Fg;

    
    Nqs = sum(nqpoints);
    
    zMat_apx = zeros(Nth, Nqs);
    zMat_apx = ANNEvalCombi(Ndof, Nr, Nth, alpha, beta, pvec, map, nqpoints, zwTable, ug_hat, zMat_apx);
    
    
figure(84);
realprt = real(zMat_apx);
mesh(xvec, yvec, realprt);    
    xlabel('Radius');
    ylabel('\theta');
    zlabel('Approximation Value');

shading interp 


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TL, Fel, xvec] = setup_LocMatTable(BdyType, left, right, alpha, beta, Nr, Nth, ...
                                        dxvec, pvec, zwTable, nqpoints, numMat, TL, Fel, vBL, vBR)


    xvec = [];
    % for each element el = 1~Nr, compute global stiffness matrix Lg and RHS fg
    Pn = 0;  % Order of each element -> Order of the last element
    for el = 1:Nr

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

        nq = nqpoints(el);
        z = zwTable(1:nq, 1, el);
        w = zwTable(1:nq, 2, el);
        xiz = (eright-eleft)/2 * z + (eleft+eright)/2;

        % Global x vector definition
        xvec = [xvec; xiz];
 
        Jn = (eright-eleft)/2;

        Nrowcol = Pn+1;
        [TL(1:Nrowcol, 1:Nrowcol, 1, el), TL(1:Nrowcol, 1:Nrowcol, 2, el), TL(1:Nrowcol, 1:Nrowcol, 3, el)] =...
            ANNElementLocal(Pn, alpha, beta, z, w, xiz, Jn);

        % Jacobian in Reparametrization variables: Depends on the function we decide when we check boundary values:
        % Note f(r,\theta) = 2 r^2 e^{r cos\theta + r sin\theta}
        % Note f(r,\theta) = [-3 r^2 + left*right] sin\theta
        Fel = ANNEvalUserLoads(left, right, alpha, beta, Pn, Jn, z, w, xiz, el, Nth, Fel);  % Real f(r,\theta) is evaluated here
        
    end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mg, Fg] = applyBdy2GMat(BdyType, left, right, Ndof, pvec, Nth, Mg, Fg, vecBL, vecBR, sigleft, sigright)


    for k = 1:Nth
        
        loc_beg = Ndof * (k - 1) + 1;
        loc_end = loc_beg + Ndof - 1;

        Mlg = Mg(loc_beg:loc_end, loc_beg:loc_end);
        Flg = Fg(loc_beg:loc_end, 1);
    
        [Mlg, Flg] = applyBdy2LMat(BdyType, left, right, Ndof, Mlg, Flg, vecBL(k, 1), vecBR(k, 1), sigleft, sigright);

        Mg(loc_beg:loc_end, loc_beg:loc_end) = Mlg;
        Fg(loc_beg:loc_end, 1) = Flg;
        
    end

return


function [Mlg, Flg] = applyBdy2LMat(BdyType, left, right, Ndof, Mlg, Flg, vBL, vBR, sigleft, sigright)

    idEnd = Ndof;

    if BdyType == 'DD'

        vecl = Lg(:, 1);      vecl(1, 1) = -1;  vecl(idEnd, 1) = 0;
        vecr = Lg(:, idEnd);  vecr(1, 1) = 0;   vecr(idEnd, 1) = -1;

        % Insert ANNecial value on Lg and Fg for boundary condition
        Mlg(1, :) = 0;       Mlg(:, 1) = 0;       Mlg(1, 1) = 1; 
        Mlg(idEnd, :) = 0;   Mlg(:, idEnd) = 0;   Mlg(idEnd, idEnd) = 1;

        Flg(1, 1) = 0; Flg(idEnd, 1) = 0;

    elseif BdyType == 'DN'

        vecl = Mlg(:, 1);       vecl(1, 1) = -1;
        vecr = zeros(Ndof, 1);  vecr(idEnd, 1) = 1;

        % Insert special value on Lg and Fg for boundary condition
        Mlg(1, :) = 0; Mlg(:, 1) = 0;   Mlg(1, 1) = 1; 

        Flg(1, 1) = 0;

    else

    end
 
    a2 = left * left;
    b2 = right * right;
    
    Flg = Flg + b2 * sigright * vecr * vBR - vecl * vBL;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zwTable, nqpoints] = setup_zwTable(qalpha, qbeta, Nr, pvec, nqpoints, zwTable)

    for ei = 1:Nr
        p_i = pvec(ei, 1);
        
        % Define z, w for gaussian quadrature
        udegree = 2 * p_i; % p+1-1 + p+1-1 = 2p
        [dgq, dgrq, dglq] = QuadratureDegree(udegree);
        
        % Gauss-Lobetto Quadrature->we take 'dglq' as a degree of corresponding Jacobi Polynomial
        
        [z, w] = JacobiGLZW(dglq, qalpha, qbeta);
        sz = size(z, 1);
        
        nqpoints(ei, 1) = sz;
        
        for j = 1:sz
            zwTable(j, 1, ei) = z(j, 1);
            zwTable(j, 2, ei) = w(j, 1);
        end
    end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coeloc = coeffglobal2local(ug_hat, map, pvec, iel)
    odr = pvec(iel, 1);
    coeloc = zeros(odr+1, 1);
    
    idxvec = map(iel, :);
    
    for i = 1 : odr+1
        k = idxvec(1, i);
        coeloc(i, 1) = ug_hat(k, 1);
    end
return
