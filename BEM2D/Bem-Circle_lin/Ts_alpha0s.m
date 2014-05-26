%   Testing the value of integral \beta_j(zi) == \int [] d/dn G(x,zi) dx

Ne = 1;
Np = 8;

us = zeros(Ne, 1); % Left circle boundary values
qs = zeros(Ne, 1); % Left circle normal values

radC = 1;
radT = 1; % This should not the same as radC

thj = 0;

thetabd = zeros(Ne, 1);
baseth = 0;

% h = 2*pi/Ne;
theta = pi/4;
h = 2*theta;
thetabd(1, 1) = - theta;

% d = h/2;

% quadrature points
qalpha = 0;
qbeta = 0;
[z, w] = JacobiGLZW(Np, qalpha, qbeta);

intval = 0; % Total sum of integral on each element

for i = 1:Ne

    thetabd(i, 1) = (i-1) * h;

    thetabd(i, 1) = - theta;
    ival = Ecalpha_Nsing(radC, thetabd(i, 1), thetabd(i, 1) + h, radT, thj, baseth, z, w);

    intval = intval + ival

end