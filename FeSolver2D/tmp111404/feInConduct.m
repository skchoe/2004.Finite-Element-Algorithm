% output 'cdtmat' is 2x2 matrix


function cdtmat = feInConduct(ei, cx, cy)

  spdim = 2;
  
  [nx, ny] = size(cx);

  szx1 = size(cx, 1);
  szx2 = size(cx, 2);
%   szy1 = size(cy, 1);
%   szy2 = size(cy, 2);
% 
%   if szx1 ~= szy1 | szx2 ~= szy2
%       TOWQUADRATURESETUNMATCHED = 0;
%   end

%     % Ex1 = \sigma(x,y) = 1;
    const = 1;
    cdtmat = subConductConst(spdim, szx1, szx2, const);

%     % Ex1.5: Nontrivial 4x4 constant tensor
%     cdtmat = subConductNtConst(spdim, szx1, szx2);

%     % Ex1.6: Nontrivial, spatially varing conductivity  (TOO BIG)
%     cdtmat = subConductFtnsp(spdim, szx1, szx2, cx, cy);

%     % Ex1.7: Nontrivial, spatially varing conductivity (diagonal)
%     cdtmat = subConductFtnspI(spdim, szx1, szx2, cx, cy);
    
%     % Ex1.8: Nontrivial, spatially varing conductivity (symmetric)
%     cdtmat = subConductFtnspII(spdim, szx1, szx2, cx, cy);

%     % Ex1.9: Nontrivial, spatially varing conductivity (symmetric+full)
%     cdtmat = subConductFtnspIII(spdim, szx1, szx2, cx, cy);

%     szx1
%     szx2
%     mat = zeros(2,2);
%     for i = 1:szx1
%         for j = 1:szx2
% 
%             mat(1, 1) = cdtmat(i, j, 1);
%             mat(1, 2) = cdtmat(i, j, 2);
%             mat(2, 1) = cdtmat(i, j, 3);
%             mat(2, 2) = cdtmat(i, j, 4);
%             xval = cx(i, j)
%             yval = cy(i, j)
%             mat
%             cholmat = chol(mat)
%         end
%     end
    
%     % Ex2 = eigenvalue = e^[x], e^[y], eigenvector = (1,1), (-1,1)
%     % Ex2 = eigen value/vectors with magnitute Mag1, Mag2
% 
%     cdtmat = subConductEigen(spdim, szx1, szx2, cx, cy);

return


function cdtmat = subConductConst(spdim, szx1, szx2, const)

    cdtmat = zeros(szx1, szx2, spdim*spdim);
    cdtmat(:,:,1) = const;
    cdtmat(:,:,4) = const;

return


function cdtmat = subConductNtConst(spdim, szx1, szx2)

    cdtmat = zeros(szx1, szx2, spdim*spdim);
    
%     nonsingmat = [2 1;3 4];
    nonsingmat = [3 50;0 1];

    matM = nonsingmat*nonsingmat';
    
%     cdtmat(:,:,1) = matM(1, 1);
%     cdtmat(:,:,2) = matM(1, 2);
%     cdtmat(:,:,3) = matM(2, 1);
%     cdtmat(:,:,4) = matM(2, 2);
    cdtmat(:,:,1) = 2;
    cdtmat(:,:,2) = 1;
    cdtmat(:,:,3) = 1;
    cdtmat(:,:,4) = 2;

return

function cdtmat = subConductFtnsp(spdim, szx1, szx2, cx, cy)

    cdtmat = zeros(szx1, szx2, spdim*spdim);

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

    cdtmat(:, :, 1) = e * e2;
    cdtmat(:, :, 2) = e24;
    cdtmat(:, :, 3) = e24;
    cdtmat(:, :, 4) = e * e22;

return

function cdtmat = subConductFtnspI(spdim, szx1, szx2, cx, cy)

    cdtmat = zeros(szx1, szx2, spdim*spdim);

    z2 = (cx.*cx + cy.*cy) * pi;
    cz2 = cos(z2)+2;

    cdtmat(:, :, 1) = cz2.*cz2;
    cdtmat(:, :, 2) = 0;
    cdtmat(:, :, 3) = 0;
    cdtmat(:, :, 4) = 1;

return


function cdtmat = subConductFtnspII(spdim, szx1, szx2, cx, cy)

    cdtmat = zeros(szx1, szx2, spdim*spdim);

    z2 = (cx.*cx + cy.*cy) * pi;
    cz21 = cos(z2);
    cz22 = cos(z2) + 2;

    cdtmat(:, :, 1) = cz21.*cz21 + cz22.*cz22;
    cdtmat(:, :, 2) = cz21;
    cdtmat(:, :, 3) = cz21;
    cdtmat(:, :, 4) = 1;

return

function cdtmat = subConductFtnspIII(spdim, szx1, szx2, cx, cy)

    cdtmat = zeros(szx1, szx2, spdim*spdim);

    z2 = (cx.*cx + cy.*cy) * pi;
    cz21 = cos(z2);
    cz22 = cz21 + 2;
    sz21 = sin(z2);
    sz22 = sz21 + 2;
    
    csz2 = cz21.*sz22;

    cdtmat(:, :, 1) = cz21.*cz21 + cz22.*cz22;
    cdtmat(:, :, 2) = csz2;
    cdtmat(:, :, 3) = csz2;
    cdtmat(:, :, 4) = sz22.*sz22;

return


function cdtmat = subConductEigen(spdim, szx1, szx2, cx, cy)

  cdtmat = zeros(szx1, szx2, spdim^2);

  for i=1:szx1
    for j=1:szx2
      cex = cx(i,j);
      cey = cy(i,j);
    
      ex = exp(cex);
      ey = exp(cey);
      cdtmat(i, j, 1) = .5*(ex+ey);
      cdtmat(i, j, 2) = .5*(ex-ey);
      cdtmat(i, j, 3) = cdtmat(i, j, 2);
      cdtmat(i, j, 4) = cdtmat(i, j, 1);
    end
  end
  

return
