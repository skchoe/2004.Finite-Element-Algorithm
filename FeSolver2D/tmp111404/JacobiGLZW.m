
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  %
  %
  % Gauss-Lobatto-Legendre
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [z,w] = JacobiGLZW(degree,alpha,beta)

  if degree == 1
    z(1) = 0.0;
    w(1) = 0.0;
  else
    one = 1.0;
    apb = alpha + beta;
    two = 2.0;
  
    z(1) = -one;
    z(degree) = one;

    if degree > 2
      z(2:degree-1) = JacobiGZeros(degree-2,alpha+one,beta+one);
    end
    w = JacobiPoly(degree-1,z,alpha,beta);
      
    fac  = (two^(apb + 1))*gamma(alpha + degree)*gamma(beta + degree);
    fac = fac/((degree-1)*gamma(degree)*gamma(alpha + beta + degree + one));
    
    w = fac./(w.*w);
    w(1) = w(1)*(beta+one);
    w(degree) = w(degree)*(alpha+one);
  end
  
  z = z';
  w = w';


