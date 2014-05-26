%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [eleft, eright, eth] = spFindIntervalByValue(left, right, dxvec, x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eleft, eright, eth] = spFindIntervalByValue(left, right, dxvec, x)

% left
% right
% dxvec

  N = size(dxvec, 1);  % Size of elements
  invs = size(x, 1);
  
  % initialize the results
  eleft = zeros(invs, 1); 
  eright = zeros(invs, 1);
  eth = zeros(invs, 1);
  
  e_th = 1;
  e_l = left;
  
  for i = 1: invs
    xe = x(i, 1);

    j1 = e_th;
    
    while (1)

      v = dxvec(j1, 1);
      e_r = e_l + v;

      if (xe >= e_l & xe < e_r) | (i == invs)
        e_th = j1;
        break;
      end

      e_l = e_r;
      j1 = j1+1;
    end

    eleft(i, 1) = e_l;
    eright(i, 1) = e_r;
    eth(i, 1) = e_th;

  end

return
