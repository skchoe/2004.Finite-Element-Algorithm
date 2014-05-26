% Error w.r.t. L^2 Norm
function l2errq = feEvalL2Error(Elm, SolVec, phiXY, Wpq, JacVec, VecXYm)

% Elm
% SolVec
% phiXY

    Nel = size(Elm, 1);
    NlocDim = size(Elm, 2);
    
    [zL, zR] = size(Wpq);

    l2err = 0;
    for ie = 1:Nel
        
        Udiffe = zeros(zL, zR);
        
        % Numerical solution
        for row = 1:NlocDim
            
            coefi = Elm(ie, row);
            phi = phiXY(:,:,row);
            sol = SolVec(coefi, 1);
            
%             if ie == 1
%                 Udiff
%             end
            Udiffe = Udiffe + (sol * phi);
%             if ie == 1
%                 sol
%                 Udiff
%             end
        end

        a11 = Udiffe(1,1);
        a21 = Udiffe(zL,1);
        a22 = Udiffe(zL, zR);

        % Getting Analytic Solution
        matx = VecXYm(:, :, ie, 1);
        maty = VecXYm(:, :, ie, 2);
        Ansol = feInSolution(matx, maty);
%         if ie == 8 | ie ==1 
%             Ansol
%             Udiffe
%         end

        Udiff1 = Udiffe - Ansol;
        Udiff2 = Udiff1 .* Udiff1;
        
        WUdiff = Wpq .* Udiff2 .* JacVec(ie, 1);

        l2err = l2err + sum(sum(WUdiff, 2));
    end
    
    
    l2errq = sqrt(l2err);
    
return
