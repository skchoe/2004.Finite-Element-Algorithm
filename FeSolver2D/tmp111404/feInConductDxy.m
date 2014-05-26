function cdtdmat = feInConductDxy(ei, cx, cy)

  spdim = 2;
    
  szx1 = size(cx, 1);
  szx2 = size(cx, 2);
  szy1 = size(cy, 1);
  szy2 = size(cy, 2);

  if szx1 ~= szy1 | szx2 ~= szy2
      TOWQUADRATURESETUNMATCHED = 0;
  end
    
%   % Ex1 = \sigma(x,y) = eye;
  const = 1;
  cdtdmat = subConductDxyConst(spdim, const, szx1, szx2, cx, cy);
    
%   % Ex1.5: Nontrivial 4x4 constant tensor
%   cdtdmat = subConductDxyNtConst(spdim, szx1, szx2, cx, cy);

% % Ex1.6:  Nontrivial, spatially varing conductivity(TOO BIG)
%     cdtdmat = subConductDxyFtnsp(spdim, szx1, szx2, cx, cy);

% % Ex1.7:  Nontrivial, spatially varing conductivity (diagonal)
%     cdtdmat = subConductDxyFtnspI(spdim, szx1, szx2, cx, cy);
  
% % Ex1.8:  Nontrivial, spatially varing conductivity (symmetric)
%     cdtdmat = subConductDxyFtnspII(spdim, szx1, szx2, cx, cy);
  
% % Ex1.9:  Nontrivial, spatially varing conductivity (symmetric+full)
%     cdtdmat = subConductDxyFtnspIII(spdim, szx1, szx2, cx, cy);


%   % Ex2 = \sigma(x,y) = x*y;
%   cdtdmat = subConductDxyEigen(spdim, szx1, szy1, cx, cy);

return


function cdtdmat = subConductDxyConst(spdim, const, szx1, szx2, cx, cy)

    cdtdmat = zeros(szx1, szx2, spdim^2);

return
  
function cdtdmat = subConductDxyNtConst(spdim, szx1, szx2, cx, cy)

    cdtdmat = zeros(szx1, szx2, spdim^2);

return

function cdtdmat = subConductDxyFtnsp(spdim, szx1, szx2, cx, cy)

    cdtdmat = zeros(szx1, szx2, spdim^2);

    z2 = cx.*cx + cy.*cy;
    z22 = z2.*z2;
    z23 = z2 + z22 + 1;
    z24 = 2*z2 + 1;

    e = exp(2)+1;
    e4 = 4*e;
    e2 = exp(2*z2);
    e22 = exp(2*z22);
    e23 = exp(z23);
    e24 = 2 * e23;
    e25 = 4*z24.*e23;

    cdtdmat(:, :, 1) = e4 * cx.*e2;
    cdtdmat(:, :, 2) = cx.*e25;
    cdtdmat(:, :, 3) = cy.*e25;
    cdtdmat(:, :, 4) = e4 * 2 * cy .* e22;
    
return

function cdtdmat = subConductDxyFtnspI(spdim, szx1, szx2, cx, cy)

    cdtdmat = zeros(szx1, szx2, spdim^2);

    z2 = (cx.*cx + cy.*cy) * pi;
    cz2 = cos(z2)+2;
    sz2 = sin(z2);

    cdtdmat(:, :, 1) = -4*pi*cx.*cz2.*sz2;
    cdtdmat(:, :, 2) = 0;
    cdtdmat(:, :, 3) = 0;
    cdtdmat(:, :, 4) = 0;
    
return

function cdtdmat = subConductDxyFtnspII(spdim, szx1, szx2, cx, cy)

    cdtdmat = zeros(szx1, szx2, spdim^2);

    z2 = (cx.*cx + cy.*cy) * pi;
    cz21 = cos(z2);
    cz22 = cos(z2) + 1;
    cz23 = cz22 + 1;
   
    sz21 = sin(z2);
    
    cdtdmat(:, :, 1) = -8*pi*cx.*sz21.*cz22;
    cdtdmat(:, :, 2) = -2*pi*cx.*sz21;
    cdtdmat(:, :, 3) = -2*pi*cy.*sz21;
    cdtdmat(:, :, 4) = 0;
    
return

function cdtdmat = subConductDxyFtnspIII(spdim, szx1, szx2, cx, cy)

    cdtdmat = zeros(szx1, szx2, spdim^2);

    z2 = (cx.*cx + cy.*cy) * pi;
    cz21 = cos(z2);
    cz22 = cos(z2) + 1;
    cz23 = cz22 + 1;
    
    c2z2 = cos(2*z2);
   
    sz21 = sin(z2);
    sz23 = sz21+2;
    csz = c2z2 - 2*sz21;
    
    
    cdtdmat(:, :, 1) = -8*pi*cx.*sz21.*cz22;
    cdtdmat(:, :, 2) = 2*pi*cx.*csz;
    cdtdmat(:, :, 3) = 2*pi*cy.*csz;
    cdtdmat(:, :, 4) = 4*pi*cy.*sz23.*cz21;
    
return

function cdtdmat = subConductDxyEigen(spdim, szx1, szx2, cx, cy)
  
  
  cdtmat = zeros(szx1, szx2, spdim^2);

  for i=1:szx1
    for j=1:szx2
      cex = cx(i,j);
      cey = cy(i,j);
    
      ex = exp(cex);
      ey = exp(cey);
      
      cdtdmat(i, j, 1) = .5*ex;
      cdtdmat(i, j, 2) = .5*ex;
      cdtdmat(i, j, 3) = -.5*ey;
      cdtdmat(i, j, 4) = .5*ey;
    end
  end
  
return


