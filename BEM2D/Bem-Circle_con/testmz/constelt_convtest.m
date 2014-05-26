% Bem with constant basis convergence test


R = 1.0;

Nbegin = 6; % Exp(2) of Nbegin
Nend = 13; %Exp(2) of Nend

xvec = zeros(Nend - Nbegin - 1, 0);
yvec = zeros(Nend - Nbegin - 1, 0);


Np = 2;    % Order of quadrature

for j = Nbegin:Nend
    
    Ne = 2^j;
    xvec(j-Nbegin+1, 1) = Ne;
    
    us = zeros(Ne, 1); % Left circle boundary values
    qs = zeros(Ne, 1); % Left circle normal values

    usol = zeros(Ne, 1); % Analytic solution

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
    
    usol;

    us = B2DCircle(Np, R, Ne, us, qs);

    diff = norm(usol - us);

    yvec(j-Nbegin+1, 1) = diff
    
end

plot(xvec, yvec);
xlabel('Num of elements');
ylabel('error to analytic solution');
title('The error graph of Boundary Element Method for constant basis');

