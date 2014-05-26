function [usol, qs] = Bm2D_CircleExact(Ne, R)

    usol = zeros(Ne, 1);

    thetabd = zeros(Ne, 1);

    h = 2*pi / Ne;
    d = h/2;

    for i = 1:Ne

        thetacnt = (i-1) * h + d;

        % Example 2 : u(x,y) = x^2 - y^2 = 0
        ctheta = cos(thetacnt);

        qs(i, 1) = 2 * R * ( 2 * ctheta * ctheta - 1 );
        usol(i, 1) = R * R * ( 2 * ctheta * ctheta - 1 );

    end

return
