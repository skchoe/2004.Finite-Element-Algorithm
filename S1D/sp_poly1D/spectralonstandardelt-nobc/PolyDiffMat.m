% compute D = differentiation matrix

  function D = PolyDiffMat(x)
  N = size(x,1);
  if N==0, D=0; return, end
  for i = 1:N,
    a(i) = 1.0;
    for j = 1:N,
      if i~=j
        a(i) = a(i)*(x(i)-x(j));
      end
    end
  end

  for i = 1:N,
    for j = 1:N,
      D(i,j) = 0.0;
      if i == j
        for k = 1:N,
          if k~=j
            D(i,j) = D(i,j) + (x(j)-x(k))^-1;
          end
        end
      else
        D(i,j) = a(i)/(a(j)*(x(i)-x(j)));
      end
    end
  end