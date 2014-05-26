% 
% The source 'f' is determined by the exact solution 'u' and Poisson equation
% - \nabula \dot \sigma\nabula u = f
% In matlab, this is implemented by
% f = - <[s11;s12;s21;s22], [uxx;uyx;uxy;uyy]> - <[s11x;s12x;s21y;s22y], [ux;uy;ux;uy]>
function value = feInSource(ei, x1, x2)

    % Prepare conductivity terms
    cdtmat = feInConduct(ei, x1, x2);
    cdtdmat = feInConductDxy(ei, x1, x2);

%     % Solution Trig1-1: Trigonometric function
% %      [ux, uy, uxx, uxy, uyx, uyy] = subSourceTrig11(x1, x2);
%     % Solution Trig1-2: Trigonometric function
%     [ux, uy, uxx, uxy, uyx, uyy] = subSourceTrig12(x1, x2);
%  
%     % Source Poly1-1: polynomial u = x^4 * y^4 with eye conductivity
%     [ux, uy, uxx, uxy, uyx, uyy] = subSourcePoly11(x1, x2);

    % Source Expo 1-3: Bounded or unbounded exponential functions
%     [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp11(x1, x2);  
%     [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp12(x1, x2);  
%    [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp13(x1, x2);  
% % % %     m1 = 0; m2 = 0; n1 = 1; n2 = 1; % Mean and std for each direction x1,x2
% % % %     [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp14(m1, m2, n1, n2, x1, x2);  

  % Solution Trig2-1: Trigonometric function to Utah Torso Data
  % [-200,200]x[-200,200]
  scaleX = 200; 
  scaleY = 200;
  [ux, uy, uxx, uxy, uyx, uyy] = subSourceTrig21(scaleX, scaleY, x1, x2);


%     % Solution Poly2-1: Polynomial to Utah Torso Data
%     % [-200,200]x[-200,200]
%   scaleX = 200; scaleY = 200;
%   [ux, uy, uxx, uxy, uyx, uyy] = subSourcePoly21(scaleX, scaleY, x1, x2);


  value = - uxx .* cdtmat(:,:,1) - ux .* cdtdmat(:,:,1)...
          - uyx .* cdtmat(:,:,2) - uy .* cdtdmat(:,:,2)...
          - uxy .* cdtmat(:,:,3) - ux .* cdtdmat(:,:,3)...
          - uyy .* cdtmat(:,:,4) - uy .* cdtdmat(:,:,4);

return


function [ux, uy, uxx, uxy, uyx, uyy] = subSourceTrig11(x1, x2)

    const11 = pi/3;
    const12 = 2*const11;
    const21 = const11^2;
    const22 = const12^2;

    Thx = x1*const11;
    Thy = x2*const12;

    cx = cos(Thx);
    sx = sin(Thx);

    cy = cos(Thy);
    sy = sin(Thy);

    ccxy = cx .* cy;
    
    ux = - const11 * sx.*cy;
    uy = - const12 * cx.*sy;

    uxx = - const21 * ccxy;
    uxy = const11 * const12 * sx .* sy;
    uyx = uxy;
    uyy = - const22 * ccxy;

return

function [ux, uy, uxx, uxy, uyx, uyy] = subSourceTrig12(x1, x2)

    const11 = pi/3;
    const12 = 2*const11;
    const21 = const11^2;
    const22 = const12^2;

    Thx = x1*const11;
    Thy = x2*const12;

    cx = cos(Thx);
    sx = sin(Thx);

    cy = cos(Thy);
    sy = sin(Thy);

    ccxy = cx .* cy;
    
    xx1 = x1.*(x1-3);
    xx2 = x2.*(x2-3);
    
    px = xx1.* sx;
    py = xx2.* sy;
    
    upx = (2*x1-3).*sx + xx1.*pi/3.*cx;
    upy = (2*x2-3).*sy + xx2 .* 2/3.*pi .* cy;

    ux = upx.*py;
    uy = px.*upy;

    uxx = ((2-const21*xx1).*sx + pi/3*(4*x1-6).* cx) .* py;
    uxy = upx.*upy;
    uyx = uxy;
    uyy = px.*((2-const22 .* xx2).*sy + 2*pi/3 *(4*x2-6).*cy);

return


function [ux, uy, uxx, uxy, uyx, uyy] = subSourcePoly11(x1, x2)

    ux = 4*x1.*x1.*x1.*x2.*x2.*x2.*x2;
    uy = 4*x1.*x1.*x1.*x1.*x2.*x2.*x2;
    uxx = 12*x1.*x1.*x2.*x2.*x2.*x2;
    uxy = 16.*x1.*x1.*x1.*x2.*x2.*x2;
    uyx = uxy;
    uyy = 12*x1.*x1.*x1.*x1.*x2.*x2;

return


function [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp11(x1, x2)
  
    Factor = exp(-x1-x2);
    ux = -Factor;
    uy = -Factor;
    uxx = Factor;
    uxy = Factor;
    uyx = Factor;
    uyy = Factor;
    
return


function [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp12(x1, x2)

    Factor = exp(x1.*x1 + x2.*x2);
    ux = 2 * x1 .* Factor;
    uy = 2 * x2 .* Factor;
    uxx = 2 * (1 + 2 * x1.*x1) .* Factor;
    uxy = 4 * (x1 .* x2) .* Factor;
    uyx = uxy;
    uyy = 2 * (1 + 2 * x2.*x2) .* Factor;
    
return

function [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp13(x1, x2)

    Factor = exp(- x1.*x1 - x2.*x2);
    ux = - 2 * x1 .* Factor;
    uy = - 2 * x2 .* Factor;
    uxx = - 2 * (1 - 2 * x1.*x1) .* Factor;
    uxy = 4 * (x1 .* x2) .* Factor;
    uyx = uxy;
    uyy = - 2 * (1 - 2 * x2.*x2) .* Factor;
    
return

function [ux, uy, uxx, uxy, uyx, uyy] = subSourceExp14(m1, m2, n1, n2, x1, x2)

    nox1 = -(x1-m1)^2 / (2*n1^2);
    nox2 = -(x2-m2)^2 / (2*n2^2);
    exy = 1/(2*pi) * 1/(n1*n2) * exp(nox1 + nox2);
    mx1 = -(x1-m1)/(n1)^2;
    mx2 = -(x2-m2)/(n2)^2;
    
    ux = mx1 .* exy;
    uy = mx2 .* exy;
    uxx = (mx1^2 - 1/n1^2).* exy;
    uxy = mx1 .* mx2 .* exy;
    uyx = uxy;
    uyy = (mx2^2 - 1/n2^2).* exy;
  
return

%-------------------------------------------------------------------
function [ux, uy, uxx, uxy, uyx, uyy] = subSourceTrig21(scaleX, scaleY, x1, x2)

    facX = 1/scaleX;
    facY = 1/scaleY;
    
    cA = facX*pi/3;
    cB = facY*2*pi/3;
    
    San = sin(cA*x1);
    Can = cos(cA*x1);
    Sbn = sin(cB*x2);
    Cbn = cos(cB*x2);
    
    
    ux = - cA .* San .* Cbn;
    uy = - cB .* Can .* Sbn;
    uxx = -cA^2 .* Can .* Cbn;
    uyx = cA*cB .* San .* Sbn;
    uxy = uyx;
    uyy = -cB^2 .* Can .* Cbn;
    
return

function [ux, uy, uxx, uxy, uyx, uyy] = subSourcePoly21(scaleX, scaleY, x1, x2)

    ux = 1/scaleX^4 * 1/scaleY^4 * 4*x1.*x1.*x1.*x2.*x2.*x2.*x2;
    uy = 1/scaleX^4 * 1/scaleY^4 * 4*x1.*x1.*x1.*x1.*x2.*x2.*x2;
    uxx = 1/scaleX^4 * 1/scaleY^4 * 12*x1.*x1.*x2.*x2.*x2.*x2;
    uxy = 1/scaleX^4 * 1/scaleY^4 * 16.*x1.*x1.*x1.*x2.*x2.*x2;
    uyx = uxy;
    uyy = 1/scaleX^4 * 1/scaleY^4 * 12*x1.*x1.*x1.*x1.*x2.*x2;
    
return
