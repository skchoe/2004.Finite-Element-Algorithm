% Ncalpha_Nsing() : Routine that evaluate the integral in terms of alpha at non-singular case
% radC
% th1, th2
% radT
% thi,
% baseth, 
% z,w

function nsa = Ealpha0_Nsing(radC, th1, th2, radT, thi, z, w)

    % Checking non singularity
    if thi == th1 | thi == th2 
        this_is_singular_case = 0;
    end


    Nzw = size(z, 1);

    xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

    fz = zeros(Nzw, 1);

    J = abs(th2-th1) / 2;

    cos_t = cos(thi);
    sin_t = sin(thi);

    bsz = zeros(Nzw, 1);    % Linear basis function term.
    
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