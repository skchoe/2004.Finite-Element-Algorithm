% Testing the value of integral 
% Result =  
Ne = 2;
Np = 4;

us = zeros(Ne, 1); % Left circle boundary values
qs = zeros(Ne, 1); % Left circle normal values

radC = 2;
radT = 0; % This should not the same as radC

thj = 3*pi/2;

thetabd = zeros(Ne, 1);

h = 2 * pi / Ne;

% quadrature points
qalpha = 0;
qbeta = 0;
[z, w] = JacobiGLZW(Np, qalpha, qbeta);

intval = 0; % Total sum of integral on each element


% thetabd(1, 1) = 0;
% ival = Ealpha0_Nsing(radC, thetabd(1,1), thetabd(1,1) + h, radT, thj, z, w)


for i = 1:Ne

    thetabd(i, 1) = (i-1) * h;

    ival = Ealpha0_Nsing(radC, thetabd(i, 1), thetabd(i, 1) + h, radT, thj, z, w);

    intval = intval + ival

end