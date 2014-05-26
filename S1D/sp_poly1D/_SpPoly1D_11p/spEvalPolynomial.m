%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function crv_img = spEvalPolynomial(alphas, Nalpha, crv_dom, Nsamp)
% alphas: array of coefficient from order 0 to Nalpha
% crv_dom: vector of domain
% Nsamp:size of crv_dom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function crv_img = spEvalPolynomial(alphas, Nalpha, crv_dom, Nsamp)

  momentMat = zeros(Nalpha, Nsamp);
  for j = 1:Nsamp
    for i = 1:Nalpha
      momentMat(i, j) = crv_dom(j, 1) ^ (i-1);
    end
  end

  img = alphas' * momentMat;
  crv_img = img';

return
