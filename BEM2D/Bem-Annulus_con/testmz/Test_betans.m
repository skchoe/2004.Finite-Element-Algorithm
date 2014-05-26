%   Testing the value of integral \beta_j(zi) == \int [] d/dn G(x,zi) dx
addpath ..

%%%%%%%%%%%%Test for integral over whole circle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ne = 50;
Np = 4;

us = zeros(Ne, 1); % Left circle boundary values
qs = zeros(Ne, 1); % Left circle normal values

radC = 2;
radT = 9; % This should not the same as radC

thj = 0;

thetabd = zeros(Ne, 1);

h = 2*pi / Ne;
% d = h/2;

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
%     ival = nonsing_beta(radC, thetabd(i, 1), thetabd(i, 1) + h, radT, thj, z, w);
% 
%     intval = intval + ival
% 
% end


radC = 2;
radT = 1;
thj = 0;
ival =Fn_beta_nonsing(radC, 0, pi/4, radT, thj, z, w)

