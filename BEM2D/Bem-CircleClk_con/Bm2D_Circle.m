function us = Bm2D_Circle(Np, R, Ne, us, qs)

    % quadrature points
    qalpha = 0;
    qbeta = 0;
    [z, w] = JacobiGLZW(Np, qalpha, qbeta);
    Nzmx = size(z, 1);

    
    % alpha, beta's
    alpha = zeros(Ne, Ne);
    beta  = zeros(Ne, Ne);
    [alpha, beta] = calcalphabeta(Np, z, w, R, Ne, alpha, beta);


    % Fs
    Fs = zeros(Ne, 1);
    Fs = calcf(Ne, alpha, beta, qs);


    % Stiffness matrix
    Lg = zeros(Ne, Ne);

    
    % Lg
    for i = 1:Ne
        for j = 1:Ne

            elt = - beta(i, j);

            if i==j
                Lg(i, j) = .5 + elt;
            else
                Lg(i, j) = elt;
            end

        end
    end

    invL = inv(Lg);
    us = invL* Fs;

return


% [alpha, beta] = calcalphabeta_RR(Np, z, w, R, Ne, alpha, beta)
function [alpha, beta] = calcalphabeta(Np, z, w, R, Ne, alpha, beta)

    Nzw = size(z, 1);

    h = 2*pi/Ne;
    d = h/2;

    for i = 1:Ne

        % angle \theta located L the center of current element.
        thi = (i-1) * h + d;

        for j = 1:Ne

            thj1 = (j-1) * h;
            thj2 = thj1 + h;

            if i==j % Integration on a singular element

                alpha(i, j) = - Alpha_sing(R, thj1, thj2, thi, z, w);

                beta(i, j)  = Beta_sing(thj1, thj2);

            else    % Integration on a non-singular element

                alpha(i, j) = - Alpha_nonsing(R, thj1, thj2, R, thi, z, w);

                beta(i, j)  = Beta_nonsing(R, thj1, thj2, R, thi, z, w);

            end

        end

    end

return



function Fs = calcf(Ne, alpha, beta, qs)

    Fs = zeros(Ne, 1);

    for i = 1:Ne

        Fs(i, 1) = dot( - alpha(i, :)', qs);

    end

return