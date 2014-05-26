%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffsine_conv.m Convergence test for spdiffsine
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spdiffsine_conv(pvec)

  N = size(pvec, 1);

  errorvec = ones(N,1);

  for i = 1:N
    errorvec(i, 1) = spdiffsine(pvec(i, 1), 0, 0, 0, 50);
  end

  semilogy(pvec, errorvec,  '.-','markersize',13);
    grid on, title('Convergence of Spectral method');
    xlabel('orders');
    ylabel('log_1_0(error)s');
    
  errorvec
  
return
