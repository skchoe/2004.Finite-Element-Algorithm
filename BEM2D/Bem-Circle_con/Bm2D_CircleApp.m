% Laplace with BEM
% r, R: inner and outer radius of annulus
% u(\theta) = r * cos (\theta)
% d/dn u(\theta) = R * sin (\theta)

% Collocation  method with locally constant basis function


R = 1.0;

k = 7;
Ne = 2^k;    % Number of element
Np = 4;    % Order of quadrature

us = zeros(Ne, 1); % Left circle boundary values
qs = zeros(Ne, 1); % Left circle normal values

usol = zeros(Ne, 1);    % Analytic solution
[usol, qs] = Bm2D_CircleExact(Ne, R);

us = Bm2D_Circle(Np, R, Ne, us, qs);

diff = norm(usol - us,2)