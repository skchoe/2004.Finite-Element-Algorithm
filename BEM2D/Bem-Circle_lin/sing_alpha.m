function nsa = sing_alpha(radC, th1, th2, thi, baseth, z, w)

    % Checking non singularity
    if thi == th1 | thi == th2 
        this_is_singular_case = 0
    end
    
    Nzw = size(z, 1);

    thj1 = abs( th1 - thi );
    thj2 = abs( th2 - thi );

    xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

    
    bsz = zeros(Nzw, 1);    % Linear basis function term.
    if baseth == 0  % \phi_0 case = 1/2 * (1+\psi);
        bsz = .5 * ( 1 + z );
    else            % \phi_1 case = 1/2 * (1-\psi);
        bsz = .5 * ( 1 - z );
    end

    
    fz1 = zeros(Nzw, 1);

    for k = 1:Nzw

        dth = abs( xtheta(k, 1) - thi ) / 2;

        if dth ~= 0
            fiz = log ( sin(dth) / dth );
        else
            fiz = 0;
        end

        fz1(k, 1) = bsz(k, 1) * fiz;

    end

    J = abs( th2 - th1 ) / 2;

    lns_term = - radC/(2*pi) * J * dot(fz, w);

    fz2 = zeros(Nzw, 1);
    
    for l = 1:Nzw
        
        fz2(l, 1) = bsz(l, 1) * log ( radC * abs(xtheta(l, 1)- thi) );
        
    end
    
    ln_term = - radC/(2*pi) * J * dot(fz2, w);

    % Sign are minus since the curve runs clockwise
    nsa = lns_term + ln_term;

return