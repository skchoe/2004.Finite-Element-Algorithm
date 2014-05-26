function nsa = Alpha_sing(radC, th1, th2, thj, z, w)

    Nzw = size(z, 1);

    thj1 = abs( th1 - thj );
    thj2 = abs( th2 - thj );

    xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

    fz1 = zeros(Nzw, 1);
    fz2 = zeros(Nzw, 1);

    for k = 1:Nzw

        dth2 = abs( xtheta(k, 1) - thj ) ;
        dth1 = dth2 / 2;

        if dth1 ~= 0
            fiz = log ( sin(dth1) ./ dth1 );
        else
            fiz = 0;
        end

        fz1(k, 1) = fiz;

    end

    J = abs( th2 - th1 ) / 2;

    lns_term = - radC/(2*pi) * J * dot(fz1, w);

    % Direct calculation for \int(log(R|x-thj|dx)
    ln_term1 = - radC/(2*pi) * ( thj2 * (log( radC*thj2 ) - 1) ...
                               + thj1 * (log( radC*thj1 ) - 1)  );

    % Sign are minus since the curve runs clockwise
    nsa = lns_term + ln_term1;

return