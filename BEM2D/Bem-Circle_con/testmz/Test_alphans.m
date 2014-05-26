% Testing the value of integral \beta_j(zi) == \int [] d/dn G(x,zi) dx
% Result =  -0.6381344945
% Ne = 1;
Np = 16;
% 
% us = zeros(Ne, 1); % Left circle boundary values
% qs = zeros(Ne, 1); % Left circle normal values
% 
% thj = 3*pi/2;
% 
% thetabd = zeros(Ne, 1);
% 
% h = pi / 2*Ne;
% % d = h/2;

addpath ../

% quadrature points
qalpha = 0;
qbeta = 0;
[z, w] = JacobiGLZW(Np, qalpha, qbeta);

radC = 10;
radT = 10;

thj = pi/32;
th1 = 0;
th2 = pi/64;

sb = Alpha_nonsing(radC, th1, th2, radT, thj, z, w)