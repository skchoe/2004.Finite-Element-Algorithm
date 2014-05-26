%%%%%%%%Evaluate f(r,theta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f is defined initially when we define the boundary values.
% Now the u(r,\theta) = -cos(\theta) -> r^2f(r,\theta) = -cos(\theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fel = ANNEvalUserLoads(left, right, alpha, beta, Pn, Jn, z, w, xiz, el, Nth, Fel)

    len = right - left;
    if len < 1e-15, exit; end
    
    Fft_f = zeros(Nth, size(xiz,1)); 

    nbgn = 0;
    nend = Nth-1;

    wNum = [nbgn : nend]';
    thetas =  wNum * 2 * pi / Nth;

    loadMat = [];
    for i = 1:size(xiz,1)
        r0 = xiz(i, 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loadF = r^2 * f(r, \theta)
%---------Ex0---------------------------------
%         loadF = ones(Nth, 1);
%---------Ex1---------------------------------        
%         loadF = cos(thetas);
%---------Ex2---------------------------------        
%         loadF = zeros(Nth, 1);
%---------Ex3---------------------------------        
        S = 2 * pi * (r0 - left) / len - pi;
        len2 = len * len;
        r02  = r0 * r0;
        pi2  = pi * pi;
        loadF = ( (r02 * 4 * pi2 / len2 + 1) * sin(S) - 2 * pi * r0 * cos(S) / len ) * cos(thetas);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%


        loadMat = [loadMat, loadF];
    end

    Fft_f = fft(loadMat);
    
    for p = 1 : Pn + 1
        Wth_p = w .* ModifiedJacobiPoly(p-1, Pn, z, alpha, beta);
        Fel(:, p, el) = 1 / Nth * Jn * Fft_f * Wth_p;
    end
    
return

