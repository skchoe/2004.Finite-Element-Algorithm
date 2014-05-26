
function alphas = spPolySCurve(prob_order)

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
  


