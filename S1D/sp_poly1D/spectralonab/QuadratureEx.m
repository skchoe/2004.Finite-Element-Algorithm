% SpDiff_SINX.m Spectral Element Solver on 1D

% Standard parametrization: = [-1, 1]

% For integral (1-x)^P*(1+x)^P on [-1, 1]
P = 20;

xdeg = [0: 1: P]';
xtimes = zeros(P+1, 1);
xterms = zeros(P+1, 1);

for p = 1:P+1
  
  % Time begins
  tic
  
  alpha = p-1;
  beta = alpha;
 
  udegree = 2*P - (alpha+beta);
  [dgq, dgrq, dglq] = QuadratureDegree(udegree);
  
  % Save # of terms in summation
  xterms(p, 1) = dglq;

  % Gauss-Lobetto Quadrature->we take 'dglq' as a degree of corresponding Jacobi Polynomial
  [z, w] = JacobiGLZW(dglq, alpha, beta);
  
  for j = 1:dglq

    zval = z(j, 1);
    wval = w(j, 1);

    fz = (1 - zval)^(P - alpha) * (1 + zval)^(P - beta);
    fmat = fmat + wval * fz;
      
  end
 
  % Time ends
  xtimes(P+1-(p-1), 1) = toc;
  
end

  subplot(2,1,1);
    plot(xdeg, xterms, '.-','markersize',13);
    xlabel('Degree of u(x)');
    grid on, title('# of Terms in Summation')

  subplot(2,1,2);
    plot(xdeg, xtimes, '.-','markersize',13);
    xlabel('Degree of u(x)');
    grid on, title('Time Elapsed')
