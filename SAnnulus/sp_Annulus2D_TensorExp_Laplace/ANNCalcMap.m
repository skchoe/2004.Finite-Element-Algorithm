%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  map = ANNCalcMap(N, pvec)
% order{1, 2, 3, ..., P, 1} version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  map = ANNCalcMap(N, pvec)

  pmax = max(pvec);

  map = zeros(N, pmax+1);
  gpinit = 1;

  gmodecount = 1;
  lmodestart = 1;
  for j = 1:N

    Pe = pvec(j, 1);
    for k = 1:Pe + 1

      if k== 1
        map(j, k) = lmodestart;
        if j~=1 gmodecount = gmodecount - 1; end
        
      elseif k== Pe+1
        map(j, k) = gmodecount;
        lmodestart = gmodecount;
        
      else
        map(j, k) = gmodecount;
        
      end

      gmodecount = gmodecount + 1;
     
    end

  end
return

