
function [alphas2, exactsolvec] = spPoissCurve(prob_order, numsampleintervals)

  alphas0 = spPolySCurve(prob_order);
  
%   % Plot of ananlytic solution
  figno = 10;
  exactsolvec = plotcurve(alphas0, prob_order, numsampleintervals, 'Plot of ananlytic solution', 'domain', 'image', figno);

  outleft0 = spEvalPolynomial(alphas0, size(alphas0,1), 0, 1)
  outright0 = spEvalPolynomial(alphas0, size(alphas0,1), 1, 1)
  
  % coefficient vector of length 1 less than order
  alphas1 = zeros(prob_order, 1);
  for i = 1:prob_order
    alphas1(i, 1) = i * alphas0(i+1, 1);
  end
%   % Plot of derivative
%   figno = 55;
%   plotcurve(alphas1, prob_order - 1, numsampleintervals+1, 'Plot of RHS force function', 'domain', 'image', figno);
  
  outleft1 = spEvalPolynomial(alphas1, size(alphas1,1), 0, 1)
  outright1 = spEvalPolynomial(alphas1, size(alphas1,1), 1, 1)
  
  % coefficient vector of length 1 less than order
  alphas2 = zeros(prob_order - 1, 1);
  for i = 1:prob_order - 1
    alphas2(i, 1) = i * (i+1) * alphas0(i+2, 1);
  end

%   % Plot of RHS
%   figno = 100;
%   plotcurve(alphas2, prob_order - 2, numsampleintervals+1, 'Plot of RHS force function', 'domain', 'image', figno);
  
return


function crv_img = plotcurve(alphas, prob_order, numsampleintervals, str_title, str_xlabel, str_ylabel, figno)

  left = 0;
  right = 1;

  length = right - left;

  unit_len = length / numsampleintervals;

  crv_dom = [left:unit_len:right]';
  crv_img = spEvalPolynomial(alphas, prob_order+1, crv_dom, numsampleintervals+1);
  
  figure(figno);
  plot(crv_dom, crv_img, '.-');
  title(str_title);
  xlabel(str_xlabel);
  ylabel(str_ylabel);
  grid on;
  
return
  
