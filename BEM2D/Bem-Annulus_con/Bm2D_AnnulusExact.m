function [usol, ui, qo, qsol] = Bm2D_AnnulusExact(Nel, R, r)

  h = 2*pi/Nel;
  d = h/2;

  thetabd = zeros(Nel, 1);
  thetacnt = zeros(Nel, 1);

  ui = zeros(Nel, 1); % Inner circle boundary values
  qsol = zeros(Nel, 1); % Inner circle normal values

  usol = zeros(Nel, 1); % Inner circle boundary values
  qo = zeros(Nel, 1); % Inner circle normal values

  for i = 1:Nel

    thetabd(i, 1) = (i-1) * h;
    thetacnt(i, 1) = (i-1) * h + d;

  end

  thetacnt
  % Example 2 : u(x,y) = x^2 - y^2 = 0
  ctheta = cos(thetacnt);

  ui = r * r * ( 2 * ctheta .* ctheta - 1 );
  qo = 2 * R * ( 2 * ctheta .* ctheta - 1 );

  usol = R * R * ( 2 * ctheta .* ctheta - 1 );
  qsol = - 2 * r * ( 2 * ctheta .* ctheta - 1 );    %Inward normal

return