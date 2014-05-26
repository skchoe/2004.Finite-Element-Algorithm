% Laplace with BEM
% r, R: inner and outer radius of annulus
% u(\theta) = r * cos (\theta)
% d/dn u(\theta) = R * sin (\theta)

% Collocation  method with locally constant basis function

R = 1;
r = 0.5;

k = 3;
Nel = 2^k; % Number of element
Np = 4;    % Order of quadrature

[usol, ui, qo, qsol] = Bm2D_AnnulusExact(Nel, R, r);

[uo, qi] = Bm2D_Annulus(Np, R, r, Nel, ui, qo);

uexsol = usol;
unumsol = uo;

u = usol./uo;

qsol;
qi;

q = qsol./qi;

diff = norm(uexsol - unumsol)/size(usol,1);
diff_max = norm(uexsol - unumsol, inf)
diff_l2 = norm(uexsol - unumsol)