
% Only visualization of solution: has nothing to do with solver, doesn't
% give exact solution to main module.
function fePlotExactSol(diV)

    len = 1/2^diV;
    [Xvec, Yvec] = meshgrid(0:len:3, 0:len:3);
    
    Nnode = 3*2^diV + 1;
    Z = zeros(Nnode, Nnode);    
    
    for i = 1:Nnode
        
        for j = 1:Nnode
            
%             if (i > 2^diV) & (i < 2*2^diV + 2) & (j > 2^diV) & (j < 2*2^diV + 2) 
%                 Z(i,j) = NaN;
%             else
                cx = len*(i-1);
                cy = len*(j-1);

                Z(i,j) = feInSolution(cx, cy);
%             end
        end
    end
    
    figure(2);
    [C,h] = contourf(Xvec, Yvec, Z, 50*(diV+1));
    
    figure(3);
    surf(Xvec, Yvec, Z);
    
return
