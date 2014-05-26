function alphas = spPolySCurveLe(prob_order)

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
  