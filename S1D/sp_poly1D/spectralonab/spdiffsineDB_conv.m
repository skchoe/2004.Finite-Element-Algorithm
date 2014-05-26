%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiffsine_conv.m Convergence test for spdiffsine
%
%   spdiffsine(left, right, samples, pvec, vdl, vdr, viz)
%
%   pvec = [1:2:15]';
%   spdiffsineDB_conv(-1, 1, 30, pvec, 0, 0, 0);
%   spdiffsineDB_conv(0, 2, 30, pvec, 1, -1, 0);
%   spdiffsineDB_conv(-50, 10, 80, pvec, -40, 20, 0);
%   spdiffsineDB_conv(-6.5, 6.5, 80, pvec, 0.5, -0.5, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left = -1;
right = 1;
samples = 90;
pvec = [1:2:15]';
vdl = 0;
vdr = 0;
viz = 0;

  N = size(pvec, 1);

  errorvec = ones(N,1);

  for i = 1:N

      [L, fvector, errorvec(i, 1)] = spdiffsineDB(left, right, samples, pvec(i, 1), vdl, vdr, viz);

  end

  semilogy(pvec, errorvec,  '.-','markersize', 13);
    grid on, title('Convergence of Spectral method');
    xlabel('orders');
    ylabel('log_1_0(error)s');
    
  errorvec
  
return
