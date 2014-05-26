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

radC = 1;
radT = 1; % This should not the same as radC

th1 = 0;
th2 = pi/8;
thj = pi/16;

ival = alpha_sing(radC, th1, th2, thj, z, w)