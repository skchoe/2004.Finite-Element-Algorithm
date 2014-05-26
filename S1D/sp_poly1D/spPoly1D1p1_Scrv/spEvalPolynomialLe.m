%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function crv_img = spEvalPolynomialLe(alphas, Nalpha, crv_dom, Nsamp)
% alphas: array of coefficient from order 0 to Nalpha
% crv_dom: vector of domain
% Nsamp:size of crv_dom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function crv_img = spEvalPolynomialLe(alphas, Nalpha, crv_dom, Nsamp)

    alpha = 0;  beta = 0; %By definition of Legendre Polynomial
    LegendreMat = zeros(Nalpha, Nsamp);
    for j = 1:Nsamp
        for i = 1:Nalpha
       
            LegendreMat(i, j) = JacobiPoly(i-1, crv_dom(j, 1), alpha, beta);
            
        end
    end

    img = alphas' * LegendreMat;
    crv_img = img';

    crv0 = crv_img(1, 1);
    crvN = crv_img(size(crv_img, 1), 1);
return
