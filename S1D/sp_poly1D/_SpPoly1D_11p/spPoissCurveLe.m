
function [crv_dom, crv_img0,crv_img1,  crv_img2] = spPoissCurveLe(prob_order, numsampleintervals)

  alphas0 = invcoeffLe(prob_order);
  
%   % Plot of ananlytic solution
  figno = 10;
  [crv_dom, crv_img0] = plotcurve(alphas0, prob_order, numsampleintervals, 'Plot of ananlytic solution', 'domain', 'image', figno);

  % Plot of derivative
  figno = 55;
  [crv_dom, crv_img1] = plotcurve1(alphas0, prob_order, numsampleintervals, 'Plot of Derivative function', 'domain', 'image', figno);

% Plot of RHS
  figno = 100;
  [crv_dom, crv_img2] = plotcurve2(alphas0, prob_order, numsampleintervals, 'Plot of RHS force function', 'domain', 'image', figno);
  
  
  figure(figno);
  plot(crv_dom, crv_img0, '.-'); hold on;
  plot(crv_dom, crv_img1, 'r.-'); hold on;  
  plot(crv_dom, crv_img2, 'g.-');
  title('S-cruve, 1, 2 Derivatives');
  legend('S-curve', 'Derivative', '2nd Derivative');
  grid on;

return


function  [crv_dom, crv_img] = plotcurve(alphas, prob_order, numsampleintervals, str_title, str_xlabel, str_ylabel, figno)

  left = 0;
  right = 1;

  length = right - left;

  unit_len = length / numsampleintervals;

  crv_dom = [left:unit_len:right]';
  crv_img = spEvalPolynomialLe(alphas, prob_order+1, crv_dom, numsampleintervals+1);
  
return
  
function  [crv_dom, crv_img] = plotcurve1(alphas, prob_order, numsampleintervals, str_title, str_xlabel, str_ylabel, figno)

  left = 0;
  right = 1;

  length = right - left;

  unit_len = length / numsampleintervals;

  crv_dom = [left:unit_len:right]';
  crv_img = spEvalPolynomialLeDeriv1(alphas, prob_order+1, crv_dom, numsampleintervals+1);

return

function  [crv_dom, crv_img] = plotcurve2(alphas, prob_order, numsampleintervals, str_title, str_xlabel, str_ylabel, figno)

  left = 0;
  right = 1;

  length = right - left;

  unit_len = length / numsampleintervals;

  crv_dom = [left:unit_len:right]';
  crv_img = spEvalPolynomialLeDeriv2(alphas, prob_order+1, crv_dom, numsampleintervals+1);

return

function alphas = invcoeffLe(prob_order)

  Nrowcol = prob_order + 1;
  
  Mat = zeros(Nrowcol, Nrowcol);

  tmp = ones(1, Nrowcol);
  hfrowcnt = 0;
  for i = 1:Nrowcol % Each row processing
        half = i/2;
        alpha = hfrowcnt;
        beta = hfrowcnt;

        if(half - floor(half) ~= 0)  %odd numbered row
            x_ = 0;
            % Order of polynomial
            for j = hfrowcnt : Nrowcol-1
                if hfrowcnt < j 
                    deg = j - hfrowcnt;
                else 
                    deg = 0;
                end
                Mat(i,j+1) = factoLe(j, hfrowcnt) * JacobiPoly(deg, x_, alpha, beta);

           end
        else
             x_  = 1;
       
            for k =  hfrowcnt : Nrowcol-1
                if hfrowcnt < k 
                    deg = k - hfrowcnt;
                else 
                    deg = 0;
                end
                Mat(i,k+1) = factoLe(k, hfrowcnt) * JacobiPoly(deg, x_, alpha, beta);

            end
            hfrowcnt = hfrowcnt + 1;
        end
  end
  
 iMat = inv(Mat);
  
  % vector of RHS when we have matrix equation.
  b = zeros(Nrowcol, 1);
  b(2,1) = 1; % Setting U(1) = 1
  
  % coefficient vector of ananlytic solution of given condition
  
  alphas = iMat * b;  

return
  
function factorial = factoLe(k, hfrowcnt)
    factorial = 1;   
    addcnt = 1;
    if hfrowcnt ~= 0
        for i = 1:hfrowcnt
            factorial = factorial * (addcnt+k)/2;
            addcnt = addcnt + 1;
        end
    end
return
