
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %
  %
  %
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function y = JacobiPolyDerivative2(degree, x, alpha, beta)

  if degree == 0 ||  degree == 1
    y = zeros(size(x));
  else
    tmp = 2.0*degree+alpha+beta;
    b1 = tmp*(1.0-x.*x);
    b2 = degree*(alpha-beta-tmp*x);
    b3 = 2.0*(degree+alpha)*(degree+beta);

    bp1 = -2*tmp.*x;
    bp2 = -degree*tmp;
    bp3 = 0;
    
    y = ( bp2.*JacobiPoly(degree,x,alpha,beta) + b2.*JacobiPolyDerivative(degree,x,alpha,beta) +...
          bp3.*JacobiPoly(degree-1,x,alpha,beta) + bp3.*JacobiPolyDerivative(degree-1,x,alpha,beta) -...
          bp1.*JacobiPolyDerivative(degree,x,alpha,beta) )./bp1;
      
  end   

