%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spStiffCurveDN.m Spectral Element Solver on 1D
%
% example commands
%  left = 0; right = 200; range of whole domain
%  dxvec = [30, 30, 40, 40, 60]'; Vector showing length of element
%  pvec = [7, 3, 5, 9, 10]'; Vector showing orders on each element
%  vdl = 0, vdr = 0: Boundary conditions are zero

% prob_order = 5
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [1, 1, 12, 1, 1]'; vdl = 0; vdr = 1; numsampleintervals = 300; %ERR(order = 5) = 0.0255
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [2, 2, 2, 2, 2]'; vdl = 0; vdr = 1; numsampleintervals = 300;  %ERR(order = 5) = 0.0019
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [3, 3, 3, 3, 3]'; vdl = 0; vdr = 1; numsampleintervals = 300;  %ERR(order = 5) = 2.4000e-004
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [4, 4, 4, 4, 4]'; vdl = 0; vdr = 1; numsampleintervals = 300;  %ERR(order = 5) = 5.6437e-006
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [5, 5, 5, 5, 5]'; vdl = 0; vdr = 1; numsampleintervals = 300;    %ERR(order = 5) = 5.7732e-015

% prob_order = 7
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [5, 5, 5, 5, 5]'; vdl = 0; vdr = 1; numsampleintervals = 300;    %ERR(order = 7) = 2.6667e-006
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [7, 7, 7, 7, 7]'; vdl = 0; vdr = 1; numsampleintervals = 300;    %ERR(order = 7) = 1.2412e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [8, 8, 8, 8, 8]'; vdl = 0; vdr = 1; numsampleintervals = 300;    %ERR(order = 7) = 1.0036e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [9, 9, 9, 9, 9]'; vdl = 0; vdr = 1; numsampleintervals = 300;    %ERR(order = 7) = 9.9920e-014
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [10, 10, 10, 10, 10]'; vdl = 0; vdr = 1; numsampleintervals = 300;    %ERR(order = 7) = 1.0636e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [11, 11, 11, 11, 11]'; vdl = 0; vdr = 1; numsampleintervals = 300;    %ERR(order = 7) = 1.0025e-013

left = 0.0;
right = 1.0;

tot_len = right - left;

%%%%%%%%%%Loop checking Error for the change of order of elements%%%%%%%%%%%
vdl = 0;
vdr = 1;

NP_ints = 5;

Perr_dom = [];
Perr_vec = [];  

NP_order_start = 6;
NP_order_end = 14;

prob_order = 10;

for ip = NP_order_start : NP_order_end

  unit_len = tot_len/NP_ints;

  dxvec = unit_len;
  pvec = ip;
  for j = 2:NP_ints
    dxvec = [dxvec;unit_len];
    pvec = [pvec; ip];
  end
  
  Psamples = 100 * ip;
  max_err_p = spStiffCurveDD(left, right, prob_order, dxvec, pvec, vdl, vdr, Psamples)

  Perr_dom = [Perr_dom, ip];
  Perr_vec = [Perr_vec, max_err_p];
  
end

figure(3);

  semilogy(Perr_dom, Perr_vec, 'o-','markersize',6, 'linewidth', 3);
    grid on;
    xlabel('order per element');
    ylabel('discrete L-infinite error');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Perr_dom = [];
% Perr_vec = [];  
% 
% NP_order_start = 1;
% NP_order_end = 11;
% 
% prob_order = 7;
% 
% for ip = NP_order_start : NP_order_end
% 
%   unit_len = tot_len/NP_ints;
% 
%   dxvec = unit_len;
%   pvec = ip;
%   for j = 2:NP_ints
%     dxvec = [dxvec;unit_len];
%     pvec = [pvec; ip];
%   end
%   
%   Psamples = 100 * ip;
%   max_err_p = spStiffCurveDN(left, right, prob_order, dxvec, pvec, vdl, vdr, Psamples)
% 
%   Perr_dom = [Perr_dom, ip];
%   Perr_vec = [Perr_vec, max_err_p];
%   
% end
% 
% figure(44);
% 
%   semilogy(Perr_dom, Perr_vec, 'o-','markersize',6, 'linewidth', 3);
%     grid on;
%     xlabel('order per element');
%     ylabel('discrete L-infinite error');
% 
