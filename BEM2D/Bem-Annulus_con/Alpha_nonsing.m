function nsa = Alpha_nonsing(radC, th1, th2, radT, thi, z, w)

    Nzw = size(z, 1);

    xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

    fz = zeros(Nzw, 1);

    J = abs(th2-th1) / 2;

    cos_t = cos(thi);
    sin_t = sin(thi);

    for k = 1:Nzw

        ck = xtheta(k, 1);
        cos_c = cos(ck);
        sin_c = sin(ck);

        expv = (radC * cos_c - radT * cos_t)^2 + (radC * sin_c - radT * sin_t)^2;
 
        fz(k, 1) = log(expv) / 2;

    end

    % Sign is minus since the parameter runs clockwise
    nsa = - radC /(2*pi) * J * dot(fz, w);

return