%   Testing the value of integral \beta_j(zi) == \int [] d/dn G(x,zi) dx

addpath '../'

% Ne = 1;
Np = 4;
% 
% us = zeros(Ne, 1); % Left circle boundary values
% qs = zeros(Ne, 1); % Left circle normal values
% 
% thetabd = zeros(Ne, 1);
% 
% % h = 2*pi/Ne;
% theta = pi/4;
% h = 2*theta;
% thetabd(1, 1) = - theta;
% 
% % d = h/2;
% 
% quadrature points
qalpha = 0;
qbeta = 0;
[z, w] = JacobiGLZW(Np, qalpha, qbeta);

th1 = 0;
th2 = pi/8;
thj = pi/16;
radC = 1;

Nzw = size(z, 1);

fz = zeros(Nzw, 1);

% Cannonical problem setting
% frac = sin(abs(z))./abs(z);
% fz = log(frac);
% 
% sb = dot(fz, w)

    thj1 = abs( th1 - thj );
    thj2 = abs( th2 - thj );

    xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

    fz1 = zeros(Nzw, 1);
    
    for k = 1:Nzw

        dth2 = abs( xtheta(k, 1) - thj ) ;
        dth1 = dth2 / 2;

        if dth1 ~= 0
            fiz = log ( sin(dth1) / dth1 );
        else
            fiz = 0;
        end

        fz1(k, 1) = fiz;

    end

    J = abs( th2 - th1 ) / 2;

    lns_term = - radC/(2*pi) * J * dot(fz1, w)

