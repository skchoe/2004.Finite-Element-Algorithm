function nsa = Fn_alpha_nonsing(radC, th1, th2, radT, thi, z, w)

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
%     nsa = - radC /(2*pi) * J * dot(fz, w);
    nsa = J * dot(fz, w);

return
% 
% function y = Fn_alpha_nonsing(x)
% 
%     Nzw = size(x, 2)
%     th1 = 0;
%     th2 = pi/8;
%     thi = pi/16;
%     
%     radC = 2;
%     radT = 1;
%     
%     cos_t = cos(thi);
%     sin_t = sin(thi);
% 
%     y = zeros(1, Nzw)
% 
%     for k = 1:Nzw
% 
%         ck = x(1, k)
%         cos_c = cos(ck);
%         sin_c = sin(ck);
% 
%         expv = (radC * cos_c - radT * cos_t)^2 + (radC * sin_c - radT * sin_t)^2;
%  
%         y(1, k) = log(expv) / 2;
% 
%     end
%     
% return
%     
%   
%     