function nsb = nonsing_beta(radC, th1, th2, radT, thj, baseth, z, w)


    % Checking non singularity
    if thi == th1 | thi == th2 
        this_is_singular_case = 0
    end
    

% Note that we have condition radC ~= radP

    Nzw = size(z, 1);

    J = abs( th2 - th1 ) / 2;

    xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

    
    bsz = zeros(Nzw, 1);    % Linear basis function term.
    if baseth == 0  % \phi_0 case = 1/2 * (1+\psi);
        bsz = .5 * ( 1 + z );
    else            % \phi_1 case = 1/2 * (1-\psi);
        bsz = .5 * ( 1 - z );
    end

    
    for k = 1:Nzw

        ck = xtheta(k, 1);

        cos_c = cos(ck);
        sin_c = sin(ck);
        
        cos_t = cos(thj);
        sin_t = sin(thj);
        
        upel = (radC * cos_c - radT * cos_t) * cos_c + (radC * sin_c - radT * sin_t) * sin_c;
        dnel = (radC * cos_c - radT * cos_t)^2 + (radC * sin_c - radT * sin_t)^2;

        fz(k, 1) = bsz(k, 1) * upel / dnel;

    end

    nsb = 1/(2*pi) * radC * J * dot(fz, w);

return