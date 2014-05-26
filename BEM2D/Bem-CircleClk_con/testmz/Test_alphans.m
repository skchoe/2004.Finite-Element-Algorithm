% Testing the value of integral \beta_j(zi) == \int [] d/dn G(x,zi) dx
% Result =  -0.6381344945
% Ne = 1;
Np = 8;
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

% intval = 0; % Total sum of integral on each element
% 
% for i = 1:Ne
% 
%     thetabd(i, 1) = (i-1) * h;
% 
%     ival = nonsing_alpha(radC, thetabd(i, 1), thetabd(i, 1) + h, radT, thj, z, w);
% 
%     intval = intval + ival
% 
% end

radC = 1;
radT = 1; % This should not the same as radC

thj = pi/2;
th1 = 0;
th2 = pi;

sb = alpha_nonsing(radC, th1, th2, radT, thj, z, w)