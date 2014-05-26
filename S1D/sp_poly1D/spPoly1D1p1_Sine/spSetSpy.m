%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L = spSetSpy(L, err, figno, viz)
%
% L: matrix to setup for Spy
% err: error threshold
% figno: figure number if this draw plot
% viz: tag that determines if we plot this or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = plotspy(L, err, figno, viz)

  Sz = size(L, 1);
  for s = 1:Sz
    for t = 1:Sz
      v = L(s,t);
      if abs(v) < err
        L(s,t) = 0;
      end
    end
  end
  
  if viz == 1
    figure(figno);
    spy(L);
  end

return