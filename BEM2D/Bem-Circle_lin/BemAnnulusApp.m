% Laplace with BEM
% r, R: inner and outer radius of annulus
% u(\theta) = r * cos (\theta)
% d/dn u(\theta) = R * sin (\theta)

% Collocation  method with locally constant basis function


R = 0.1;

Ne = 64;    % Number of element
Np = 2;    % Order of quadrature

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

diff = norm(usol - us)

