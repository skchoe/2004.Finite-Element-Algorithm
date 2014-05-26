function crv_img = spEvalPolynomialLeDeriv2(alphas, Nalpha, crv_dom, Nsamp)

    alpha = 2;  beta = 2; %By definition of 2nd Deriv. of Legendre Polynomial
    alphas2 = alphas(3:size(alphas, 1), 1);
    
    LegendreMat = zeros(Nalpha-2, Nsamp);
    
    for j = 1:Nsamp
        for i = 1:Nalpha-2
 
            LegendreMat(i, j) = 0.25 * (i+2) * (i+3) * JacobiPoly(i-1, crv_dom(j, 1), alpha, beta);
            
        end
    end

    img = alphas2' * LegendreMat;
    crv_img = img';

    drv2_0 = crv_img(1, 1);
    drv2_N = crv_img(Nsamp, 1);
    
return
