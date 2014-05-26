
function [alphas2, exactsolvec] = spPoissCurve(prob_order, numsampleintervals)

  alphas0 = invcoeff(prob_order);
  
%   % Plot of ananlytic solution
  figno = 10;
  exactsolvec = plotcurve(alphas0, prob_order, numsampleintervals, 'Plot of ananlytic solution', 'domain', 'image', figno);

  
%   % Plot of derivative
%   figno = 55;
%   plotcurve(alphas1, prob_order - 1, numsampleintervals+1, 'Plot of RHS force function', 'domain', 'image', figno);
  
  
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
  

function alphas = invcoeff(prob_order)

  Nrowcol = prob_order + 1;
  
  Mat = zeros(Nrowcol, Nrowcol);
%   iMat = inv(Mat);

  tmp = ones(1, Nrowcol);
  rowcnt = 1;
  for i = 1:Nrowcol
    half = i/2;
    if(half - floor(half) ~= 0)  %odd numbered row
      col = (i+1)/2;
      Mat(i, col) = facto(col-1);
    else
      if i == 2
        Mat(i, :) = tmp(1, :);
      else
        for j = 1:Nrowcol
          Mat(i, j) = tmp(1, j) * (j-rowcnt);
        end
        tmp(1, :) = Mat(i, :);
        rowcnt = rowcnt + 1;
      end
    end
  end

  miniM = min(min(Mat));
  maxiM = max(max(Mat));
  
  iMat = inv(Mat);
  
  miniIM = min(min(iMat));
  maxiIM = max(max(iMat));
  
  % vector of RHS when we have matrix equation.
  b = zeros(Nrowcol, 1);
  b(2,1) = 1; % Setting U(1) = 1
  
  % coefficient vector of ananlytic solution of given condition
  
  alphas = iMat * b;  

return
  
function factorial = facto(n)
  factorial = 1;
  for i = 1:n
    factorial = factorial * (n-i+1);
  end
return
  


