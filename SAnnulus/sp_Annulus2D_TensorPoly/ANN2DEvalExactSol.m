%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xvec, yvec_sol] = ANN1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xvec, yvec, zMat_sol] = ANN2DEvalExactSol(left, right, xvec, yvec)

    len = right - left;
    if len == 1e-15, exit; end
    
    Nqr = size(xvec, 1);
    Nth = size(yvec, 1);
    
    zMat_sol = zeros(Nth, Nqr);
    
    for i = 1:Nqr
        % Ex1 : zMat_sol(:, i) = -cos(yvec);                    
        % Ex2 : zMat_sol(:, i) = - xvec(i, 1) * cos(yvec);
        S = 2 * pi * (xvec(i, 1) - left) / len - pi;
        zMat_sol(:, i) = sin(S) * cos(yvec);
    end

%   %--------------------------------------------------------------  
    figno = 21;
    figure(figno);
  

    % Analytic solution considering boundary condition
%     mesh(zMat_sol)
%     mesh(xvec, yvec, zMat_sol);    
    surf(xvec, yvec, zMat_sol);    
        xlabel('Radius');
        ylabel('\theta');
        zlabel('Exact Solution');

    shading interp


return
    
% function [xvec, yvec, zMat_sol] = ANN2DEvalExactSol(BdyType, left, right, xvec, yvec, vBL, vBR);
% 
%     nx = size(xvec, 1);
%     ny = size(yvec, 1);
% 
%     zMat_sol = zeros(ny, nx);
%     
%   
%     if  BdyType == 'DN'
%         uls = ANNEvalUserSolution(left, right, left, ny, 0)';
%         urs = ANNEvalUserSolution(left, right, right, ny, 1)';
%         
%         C = vBR - urs;
%         D = vBL - uls - C * left;
% 
%     
%     elseif BdyType == 'DD'
%         uls = ANNEvalUserSolution(left, right, left, ny, 0)';
%         urs = ANNEvalUserSolution(left, right, right, ny, 0)';
% 
%         d1 = - vBR + vBL;
%         e1 = urs - uls;
%         dlen = right - left;
%     
%         C = (d1 - e1) / dlen;
%         D = - vBR - urs - C * right;
%                                                      
%     else
%     end
% 
%     uvec = ANNEvalUserSolution(left, right, xvec, ny, 0);
% 
%     for i = 1:ny
%         zMat_sol(i, :) = uvec(:, i)' + C(i, 1) * xvec' + D(i, 1)*ones(nx, 1)';
%     end
%     
% %     zMat_sol
% %   %--------------------------------------------------------------  
%     figno = 21;
%     figure(figno);
%   
%     % Analytic solution considering boundary condition
%     mesh(xvec, yvec, zMat_sol);
%     shading interp 
% 
% 
% return
