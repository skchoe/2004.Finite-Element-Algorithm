

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

function zMat_apx = ANNEvalCombi(Ndof, Nr, Nth, alpha, beta, pvec, map, nqpoints, zwTable, ug_hat, zMat_apx);


    Nqs = sum(nqpoints);

    % Reorder modes for its angular modes and quadrature orders.
    ie = 1;
    for ith = 1:Nth
        for j = 1:Ndof
            U(ith, j) = ug_hat(ie, 1);
            ie = ie+1;
        end
    end

    uHatMat = [];
    for kf = 1:Ndof
        w = U(:, kf);
        uHatMat = [uHatMat ifft(w)];
    end
    
    % Basis Value Table Decleration
    bsTbMat = zeros(Nqs, Ndof);
    bsTbMat = setup_basisTables(alpha, beta, Nr, nqpoints, pvec, map, zwTable, bsTbMat);


    basis = bsTbMat';
    zMat_apx = uHatMat * bsTbMat' * Nth;

    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bsTbMat = setup_basisTables(alpha, beta, Nr, nqpoints, pvec, map, zwTable, bsTbMat)

    rowinfo = 0;
    for ei = 1:Nr
        p_i = pvec(ei, 1);
        z_i = nqpoints(ei, 1);

        startId = map(ei, 1);
        colId = startId;
        
        for ip = 0:p_i
            for iz = 1:z_i
                z = zwTable(iz, 1, ei);
                bsTbMat(rowinfo + iz, colId) = ModifiedJacobiPoly(ip, p_i, z, alpha, beta);
            end
            colId = colId + 1;
        end
        
        rowinfo = rowinfo + z_i;
    end

return
