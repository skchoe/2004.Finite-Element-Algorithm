%--------------------------------------------------------------------------------------
function [w, etime] = spectraldiffbymatrix(v, N, h, bt)

  if bt ~= 0
      tic;
  end
  
  column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
%   colsize = size(column)
%   vsize = size(v)
  D = toeplitz(column,column([1 N:-1:2]));
  w1 = D*v';
  w = w1';
  
  if bt ~= 0
    etime = toc;
  else 
    etime = 0;
  end
  
return
