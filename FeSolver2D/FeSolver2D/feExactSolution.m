% Ndm : 1: Node index, 2: x-Coordinate 3: y-coordinate

function [solIn, solBd] = feExactSolution(Ndm, BdyN, IntNd)

    Nin = size(IntNd, 1);
    solIn = zeros(Nin, 1);
    
    for i=1:Nin
        cidx = IntNd(i, 1);
        cx = Ndm(cidx, 2);
        cy = Ndm(cidx, 3);

        solex = feInSolution(cx, cy);
        solIn(i, 1) = solex;
    end
    
    
    NbdN = size(BdyN, 1);
    solBd = zeros(NbdN, 1);

    if NbdN > 1
        for j = 1:NbdN
            cidx = BdyN(j, 1);
            cx = Ndm(cidx, 2);
            cy = Ndm(cidx, 3);
        
            solex = feInSolution(cx, cy);
            solBd(j, 1) = solex;
        end
    end

return
