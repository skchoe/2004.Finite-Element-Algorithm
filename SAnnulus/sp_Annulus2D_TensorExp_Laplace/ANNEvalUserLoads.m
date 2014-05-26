%%%%%%%%Evaluate f(r,theta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f is defined initially when we define the boundary values.
% Now the u(r,\theta) = -cos(\theta) -> r^2f(r,\theta) = -cos(\theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fel = ANNEvalUserLoads(left, right, alpha, beta, Pn, Jn, z, w, xiz, el, Nth, Fel)

    len = right - left;
    if len < 1e-15, exit; end
    
    Fft_f = zeros(Nth, size(xiz,1)); 

    nbgn = 0;
    nend = Nth-1;

    wNum = [nbgn : nend]';
    thetas =  wNum * 2 * pi / Nth;

    loadMat = [];
    for i = 1:size(xiz,1)
        r0 = xiz(i, 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loadF = r^2 * f(r, \theta)
%-------Ex0---------------------------------
%         loadF = ones(Nth, 1);
%-------Ex1---------------------------------        
%          loadF = cos(thetas);
%-------Ex2---------------------------------        
%         loadF = zeros(Nth, 1);
%-------Ex3---------------------------------        
%        S = 2 * pi * (r0 - left) / len - pi;
        
%        ten0 = ANNEvalUserTensor(1, r0);
%        ten1 = ANNEvalUserTensorDerivative(1, r0);
%        
%        len2 = len * len;
%        pi2 = pi * pi;
%        
%        u_r = 2 * pi * cos(S) * cos(thetas) / len;
%        u_rr = - 4 * pi2 * sin(S) * cos(thetas) / len2;
%        u_thth = - sin(S) * cos(thetas);
%        
%        loadF = - r0 * r0 * ( ten1 * u_r + ten0 * u_rr ) - r0 * ten0 * u_r - ten0 * u_thth;
%-------Ex4---------------------------------        
%         
%         prob_order = 9;
%         alphas = spPolySCurveLe(prob_order);
%         Nalpha = size(alphas, 1);
%         
%         exp_coef = exp(-(r0-1)*(r0-1));
%         
%         deriv1_coef = 2*r0*r0*(r0-1)-r0;
%         deriv2_coef = -r0*r0;
%         
%         deriv0_term = spEvalPolynomialLe(alphas, Nalpha, r0, 1);
%         deriv1_term = spEvalPolynomialLeDeriv1(alphas, Nalpha, r0, 1);
%         deriv2_term = spEvalPolynomialLeDeriv2(alphas, Nalpha, r0, 1);
%         
%         
%         loadF = exp_coef * ( deriv1_coef * deriv1_term + deriv2_coef * deriv2_term + deriv0_term ) * cos(thetas);
%-------Ex5 : Laplace Eqn ---------------------------------        
%         
        loadF = zeros(Nth, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%


        loadMat = [loadMat, loadF];
    end

    Fft_f = fft(loadMat);
    
    for p = 1 : Pn + 1
        Wth_p = w .* ModifiedJacobiPoly(p-1, Pn, z, alpha, beta);
        Fel(:, p, el) = 1 / Nth * Jn * Fft_f * Wth_p;
    end
    
return

