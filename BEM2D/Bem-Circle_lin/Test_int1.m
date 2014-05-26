% Nonsing integral of log(sin(x)/x).

Ne = 1;
Np = 8;

us = zeros(Ne, 1); % Left circle boundary values
qs = zeros(Ne, 1); % Left circle normal values

% quadrature points
qalpha = 0;
qbeta = 0;
[z, w] = JacobiGLZW(Np, qalpha, qbeta);

left =pi/4;
right = 3*pi/4;

xiz = (right-left)/2 * z + (right+left)/2

J = abs(right-left)/2;

%Quadrature formula

Nzw = size(z, 1);

fiz = zeros(Nzw, 1);

for i = 1:Nzw
    x = xiz(i, 1);
    fiz(i, 1) = log( sin(x)/x );
end

intval = J * dot(fiz, w)
