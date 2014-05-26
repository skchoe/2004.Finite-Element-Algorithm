% JacobiPolyDerivativeN.m : Given degree, order = N, alpha, beta, indep.var.x,
%                           Calculate Y = d^order/dx^order * p(x)
function y = JacobiPolyDerivativeN(degree, x, alpha, beta, order)

  if order == 0
      Y = JacobiPoly(degree, x, alpha, beta);
  else
    if degree <= order - 1
      y = zeros(size(x));
    else
      coeff = 1;
      for i = 1:order
        tmp = 0.5 * (alpha + i-1 + beta + i-1 + degree - i + 1 + 1);
        coeff = coeff * tmp;
      end
      y = coeff * JacobiPoly(degree - order, x, alpha + order, beta + order);
    end
  end

