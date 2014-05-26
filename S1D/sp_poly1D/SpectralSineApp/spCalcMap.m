%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  map = spCalcMap(N, pvec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  map = spCalcMap(N, pvec)

  pmax = max(pvec);

  map = zeros(N, pmax+1);
  gpinit = 1;

  gmodecount = 1;
  lmodestart = 1;
  for j = 1:N

    for k = 1:pvec(j, 1) + 1

      if k== 1
        map(j, k) = lmodestart;
        if j~=1 gmodecount = gmodecount - 1; end
      elseif k== 2
        map(j, k) = gmodecount;
        lmodestart = gmodecount;
      else
        map(j, k) = gmodecount;
      end

      gmodecount = gmodecount + 1;
     
    end

  end
return

