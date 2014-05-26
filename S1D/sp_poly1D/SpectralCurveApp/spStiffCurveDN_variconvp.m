%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spStiffCurveDN.m Spectral Element Solver on 1D
%
% example commands
%  left = 0; right = 200; range of whole domain
%  dxvec = [30, 30, 40, 40, 60]'; Vector showing length of element
%  pvec = [7, 3, 5, 9, 10]'; Vector showing orders on each element
%  vdl = 0, vnr = 0: Boundary conditions are zero

% prob_order = 5
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [1, 1, 12, 1, 1]'; vdl = 0; vnr = 0; numsampleintervals = 300; %ERR(order = 5) = 0.0255
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [2, 2, 2, 2, 2]'; vdl = 0; vnr = 0; numsampleintervals = 300;  %ERR(order = 5) = 0.0019
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [3, 3, 3, 3, 3]'; vdl = 0; vnr = 0; numsampleintervals = 300;  %ERR(order = 5) = 2.4000e-004
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [4, 4, 4, 4, 4]'; vdl = 0; vnr = 0; numsampleintervals = 300;  %ERR(order = 5) = 5.6437e-006
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [5, 5, 5, 5, 5]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 5) = 5.7732e-015

% prob_order = 7
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [5, 5, 5, 5, 5]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 2.6667e-006
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [7, 7, 7, 7, 7]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.2412e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [8, 8, 8, 8, 8]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.0036e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [9, 9, 9, 9, 9]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 9.9920e-014
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [10, 10, 10, 10, 10]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.0636e-013
% left = 0; right = 1; dxvec = [.2, .2, .2, .2, .2]'; pvec = [11, 11, 11, 11, 11]'; vdl = 0; vnr = 0; numsampleintervals = 300;    %ERR(order = 7) = 1.0025e-013

left = 0.0;
right = 1.0;

tot_len = right - left;

vdl = 0;
vnr = 0;

NH_ints = 5;
unit_len = tot_len/NH_ints;

dxvec = unit_len;
for j = 2:NH_ints
  dxvec = [dxvec;unit_len];
end

%%%%%%%%%%Loop checking Error for the change of order of elements%%%%%%%%%%%

prob_order = 5;
NH_order_start = 1;
NH_order_end = (prob_order+1)/2;
N_curves = NH_order_start - NH_order_end + 1;

NP_order_start = 1;
NP_order_end = 5;
N_orders = NP_order_start - NP_order_end + 1;

Perr_vec = zeros(N_orders, N_curves);


%NH_order_end = 1;
for ih = NH_order_start : NH_order_end
%ih = 1, 2, 3->2k-1:1, 3, 5
  
  ihnew = 2*ih - 1;
  Perr_dom = [];
  
  for ip = NP_order_start : NP_order_end
    
    pvec = [ihnew;ihnew;ip;ihnew;ihnew];
  
    Psamples = 100 * ihnew;
    max_err_p = spStiffCurveDN(left, right, prob_order, dxvec, pvec, vdl, vnr, Psamples);

    Perr_dom = [Perr_dom; ip];
    Perr_vec(ip, ih) = max_err_p;
    
  end
  
end
  
figure(33);
  
  semilogy(Perr_dom, Perr_vec(:, 1), 'o-',  'markersize',6, 'linewidth', 3); hold on;
  semilogy(Perr_dom, Perr_vec(:, 2), 's-.', 'markersize',6, 'linewidth', 3); hold on;
  semilogy(Perr_dom, Perr_vec(:, 3), 'v--', 'markersize',6, 'linewidth', 3); hold off;
    grid on;
    xlabel('order per element');
    ylabel('discrete L-infinite error');
    legend('order[1, 1, p, 1, 1]','order[3, 3, p, 3, 3]','order[5, 5, p, 5, 5]',3);

hold off;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

prob_order = 7;
NH_order_start = 1;
NH_order_end = (prob_order+1)/2;
N_curves = NH_order_start - NH_order_end + 1;

NP_order_start = 1;
NP_order_end = 7;
N_orders = NP_order_start - NP_order_end + 1;

Perr_vec = zeros(N_orders, N_curves);


for ih = NH_order_start : NH_order_end
%ih = 1, 2, 3, 4 -> 2k-1:1, 3, 5, 7
  
  ihnew = 2*ih - 1;
  Perr_dom = [];
  
  for ip = NP_order_start : NP_order_end
    
    pvec = [ihnew;ihnew;ip;ihnew;ihnew];
  
    Psamples = 100 * ihnew;
    max_err_p = spStiffCurveDN(left, right, prob_order, dxvec, pvec, vdl, vnr, Psamples);

    Perr_dom = [Perr_dom; ip];
    Perr_vec(ip, ih) = max_err_p;
    
  end
  
end
  
figure(44);

  semilogy(Perr_dom, Perr_vec(:, 1), 'o-',  'markersize',6, 'linewidth', 3); hold on;
  semilogy(Perr_dom, Perr_vec(:, 2), 's-.', 'markersize',6, 'linewidth', 3); hold on;
  semilogy(Perr_dom, Perr_vec(:, 3), 'v--', 'markersize',6, 'linewidth', 3); hold on;
  semilogy(Perr_dom, Perr_vec(:, 4), '*:', 'markersize',6, 'linewidth', 3); hold off;
    grid on;
    xlabel('order per element');
    ylabel('discrete L-infinite error');
    legend('order[1, 1, p, 1, 1]','order[3, 3, p, 3, 3]','order[5, 5, p, 5, 5]','order[7, 7, p, 7, 7]',3);


hold off;
