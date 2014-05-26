%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffsine_conv.m Convergence test for spdiffsine
%
%   spdiffsine(left, right, samples, pvec, vdl, vnr, viz)
%
%   pvec = [1:2:15]';
%   spdiffsineNB_conv(-0.5, 0.3, 30, pvec, -1, 100, 1);
%   spdiffsineNB_conv(-1, 1, 30, pvec, 0, 0, 1);
%   spdiffsineNB_conv(0, 2, 30, pvec, 1, -1, 1);
%   spdiffsineNB_conv(-5, 10.5, 80, pvec, -40, 10, 1);
%   spdiffsineNB_conv(-6.5, 6.5, 80, pvec, 0.5, -0.5, 1);
%   spdiffsineNB_conv(100.5, 101.4, 40, pvec, 50.5, -1.5, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spdiffsineNB_conv(left, right, samples, pvec, vdl, vnr, viz)

  N = size(pvec, 1);

  errorvec = ones(N,1);

  for i = 1:N
    [L, fvector, errorvec(i, 1)] = spdiffsineNB(left, right, samples, pvec(i, 1), vdl, vnr, viz);
  end
  
  
figure(3);
  semilogy(pvec, errorvec,  '.-','markersize',13);
    grid on, title('Convergence of Spectral method');
    xlabel('orders');
    ylabel('log_1_0(error)s');
    
  errorvec
  
return
