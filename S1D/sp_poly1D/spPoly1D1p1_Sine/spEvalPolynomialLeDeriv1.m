function crv_img = spEvalPolynomialLeDeriv1(alphas, Nalpha, crv_dom, Nsamp)

    alpha = 1;  beta = 1; %By definition of 2nd Deriv. of Legendre Polynomial
    alphas1 = alphas(2:size(alphas, 1), 1);
    
    LegendreMat = zeros(Nalpha-1, Nsamp);
    
    for j = 1:Nsamp
        for i = 1:Nalpha-1
 
            LegendreMat(i, j) = 0.5 * (i+1) * JacobiPoly(i-1, crv_dom(j, 1), alpha, beta);
            
        end
    end

    img = alphas1' * LegendreMat;
    crv_img = img';

    drv1_0 = crv_img(1, 1);
    drv1_N = crv_img(Nsamp, 1);
    
return
