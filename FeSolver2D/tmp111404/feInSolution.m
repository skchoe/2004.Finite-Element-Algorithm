function value = feInSolution(x1, x2)

%     % Solution Trig 1-1: Trigonometric function
% %      value = subSolutionTrig11(x1, x2);
%     value = subSolutionTrig12(x1, x2);  %Consine * x * (x-3) so that this satifies general Neumann condition dot(n, \sig\del u) = 0;
%     
%     % Solution Poly 1-2: Polynomial
%     value = subSolutionPoly11(x1, x2);

%     Solution Expo 1-3: Bounded or unbounded exponential functions
%     value = subSolutionExp11(x1, x2);  
%     value = subSolutionExp12(x1, x2);  
%     value = subSolutionExp13(x1, x2);  
% % % %     m1 = 0; m2 = 0; n1 = 1; n2 = 1; % Mean and std for each direction x1,x2
% % % %     value = subSolutionExp14(m1, m2, n1, n2, x1, x2);  


%     Solution Trig 2-1: Trigonometric function to Utah Torso Data
%     [-200,200]x[-200,200]
    scaleX = 200;
    scaleY = 200;
    value = subSolutionTrig21(scaleX, scaleY, x1, x2);

%     % Solution Poly 2-1: Polynomial to Utah Torso Data
%     % [-200,200]x[-200,200]
%     scaleX = 200;
%     scaleY = 200;
%     value = subSolutionPoly21(scaleX, scaleY, x1, x2);

return

function value = subSolutionTrig11(x1, x2)
    
    value = cos(x1*pi/3) .* cos(x2*2*pi/3);

return

function value = subSolutionTrig12(x1, x2)
    
    value = x1.*(x1-3).*sin(x1*pi/3) .* x2.*(x2-3).*sin(x2*2*pi/3);

return

function value = subSolutionPoly11(x1, x2)
    
    value = x1 .* x1 .* x1 .* x1 .* x2 .* x2 .* x2 .* x2;
    
return

function value = subSolutionExp11(x1, x2)

    value = exp(-x1-x2);

return


function value = subSolutionExp12(x1, x2)

    value = exp(x1.*x1 + x2.*x2);
    
return


function value = subSolutionExp13(x1, x2)

    value = exp(- x1.*x1 - x2.*x2);
    
return

function value = subSolutionExp14(m1, m2, n1, n2, x1, x2)

    nox1 = -(x1-m1)^2 / (2*n1^2);
    nox2 = -(x2-m2)^2 / (2*n2^2);
    value = 1/(2*pi) * 1/(n1*n2) * exp(nox1 + nox2);

return

function value = subSolutionTrig21(scaleX, scaleY, x1, x2)
    
    facX = 1/scaleX;
    facY = 1/scaleY;
    x1 = facX * x1;
    x2 = facY * x2;
    value = cos(x1*pi/3) .* cos(x2*2*pi/3);

return

function value = subSolutionPoly21(scaleX, scaleY, x1, x2)

    x1 = 1/scaleX * x1;
    x2 = 1/scaleY * x2;
    value = x1 .* x1 .* x1 .* x1 .* x2 .* x2 .* x2 .* x2;
    
return

