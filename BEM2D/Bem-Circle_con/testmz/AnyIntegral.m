Np = 8;
addpath ../

% quadrature points
qalpha = 0;
qbeta = 0;
[z, w] = JacobiGLZW(Np, qalpha, qbeta);

radC = 1;
radT = .5; % This should not the same as radC

thj = pi/16;
th1 = 0;
th2 = pi/8;

J = abs(th2-th1) / 2;

xtheta = (th2 - th1) / 2 * z + (th2 + th1) / 2;

Nz = size(z,1);

fz1 = zeros(Nz,1);
fz2 = zeros(Nz,1);
for k = 1:Nz
  x = xtheta(k,1);
%   th = thj - x;
%   fz1(k, 1) = log(radC^2 + radT^2 - 2*radC*radT*th);
  fz1(k, 1) = log(cos(x));
  fz2(k, 1) = x * tan(x);
end

    % Sign is minus since the parameter runs clockwise
nsa1 = - radC /(4*pi) * J * dot(fz1, w)

nsa2 = th2*log(cos(th2)) - th1*log(cos(th1)) - J*dot(fz2,w)
