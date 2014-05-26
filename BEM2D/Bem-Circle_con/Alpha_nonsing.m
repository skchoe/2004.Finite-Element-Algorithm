% This compute the following integral
%
% \alpha(x) = - radC/(2*Pi) * int(0.5*log(|z-x|^2), x=th1...th2);
%
%

function nsa = Alpha_nonsing(radC, th1, th2, radT, thi, z, w)

%   if radC ~= radT

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


    
% Only when two radii are same
% -> - R/(2*pi) * int(Rl|x-x0|) x=th1...th2
% %   else
% 
%     R = radC;
% 
%     d2 = abs(th2 - thi);
%     d1 = abs(thi - th1);
%
%     log2 = log(R*d2);
%     log1 = log(R*d1);
%
%     nsa = -R/(2*pi) * ( d2 * (log2 - 1) + d1 * (log1 - 1) )
%
%     end

return