function sb = sing_beta(radC, th1, th2, thi, baseth, z, w)

    % Checking non singularity
    if thi == th1 | thi == th2 
        this_is_singular_case = 0
    end
    
    J = abs( th1 - th2 )/2;
    
    Nzw = size(z, 1);
    
    bsz = zeros(Nzw, 1);    % Linear basis function term.
    if baseth == 0  % \phi_0 case = 1/2 * (1+\psi);
        bsz = .5 * ( 1 + z );
    else            % \phi_1 case = 1/2 * (1-\psi);
        bsz = .5 * ( 1 - z );
    end

    sb = 1/(4*pi) * J * dot(bsz, w);
    
return