%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffpoly_conv.m Convergence test for spdiffpoly
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spdiffpoly_conv(pvec)

  N = size(pvec, 1);

  errorvec = ones(N,1);

  for i = 1:N
    errorvec(i, 1) = spdiffpoly(pvec(i, 1), 0);
  end

  plot(pvec, errorvec,  '.-','markersize',13);
    grid on, title('Convergence of Spectral method');
    xlabel('orders');
    ylabel('errors');
  errorvec
  

  
return
