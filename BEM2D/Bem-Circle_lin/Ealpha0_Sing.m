function nsa = sing_alpha(radC, th1, th2, thi, z, w)

    Nzw = size(z, 1);

    thj1 = abs( th1 - thi );
    thj2 = abs( th2 - thi );

    xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

    fz = zeros(Nzw, 1);

    for k = 1:Nzw

        dth = abs( xtheta(k, 1) - thi ) / 2;

        if dth ~= 0
            fiz = log ( sin(dth) / dth );
        else
            fiz = 0;
        end

        fz(k, 1) = fiz;

    end

    J = abs( th2 - th1 ) / 2;

    lns_term = - radC/(2*pi) * J * dot(fz, w);

%     logand1 = radC * thj2
%     logand2 = radC * thj1
%     
    ln_term = - radC/(2*pi) * ( radC * thj2 * (log( radC*thj2 ) - 1) + radC * thj1 * (log( radC*thj1 ) - 1));

    % Sign are minus since the curve runs clockwise
    nsa = lns_term + ln_term;

return