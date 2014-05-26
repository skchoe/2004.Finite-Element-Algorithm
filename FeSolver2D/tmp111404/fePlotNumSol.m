% function fePlotNumSol(figno, plot_title, Ndm, IntNd, BdyN, BdyD, ...
%		      SolIn, SolBd, BdyDP, Elm) 
%
%
% inputs:
%   figno- the number of the figure for plotting
%   bsOdr- 
%   Elm- the connectivities of the elements
%   Ndm- the geometric locations of the points



function fePlotNumSol(figno, plot_title, bsOdr, Ndm, IntNd, BdyN, BdyD, ...
		      SolIn, SolBd, BdyDP, Elm) 

  Nnd = size(Ndm, 1);
  Xvec = zeros(Nnd, 1);
  Yvec = zeros(Nnd, 1);
  Sovec = zeros(Nnd, 1);

 %sizenddm = size(Ndm)
 
  % Contains All info(Total Nodes, Total BdyD/N)
  Nin = size(IntNd, 1); 
  for i = 1:Nin
      %arg = IntNd(i, 1)
    Xvec(i, 1) = Ndm(IntNd(i, 1), 2);
    Yvec(i, 1) = Ndm(IntNd(i, 1), 3);
    Sovec(i, 1) = SolIn(i, 1);  
  end
  
  NbdN = size(BdyN, 1);
  for i = 1:NbdN
    Xvec(Nin + i, 1) = Ndm(BdyN(i, 1), 2);
    Yvec(Nin + i, 1) = Ndm(BdyN(i, 1), 3);
    Sovec(Nin+ i, 1) = SolBd(i, 1);
  end
  
  NbdD = size(BdyD, 1);
  for i = 1:NbdD
    Xvec(Nin+NbdN+i, 1) = Ndm(BdyD(i, 1), 2);
    Yvec(Nin+NbdN+i, 1) = Ndm(BdyD(i, 1), 3);
    Sovec(Nin+NbdN+i, 1) = BdyDP(i, 1);
  end

  
  % Compose Element based on approximation order bsOdr
  if bsOdr == 'L'
    NgeoElm = Elm(:,1:3);
  elseif bsOdr == 'Q'
    NgeoElm = subCompNewEltsQuadratic(Elm);
  elseif bsOdr == 'C'
    NgeoElm = subCompNewEltsCubic(Elm);
  end
  
  figure(figno);
%   trimesh(Elm, Xvec, Yvec);
  trisurf(NgeoElm, Xvec, Yvec, Sovec);
  title(plot_title);
    
return

function NgeoElm = subCompNewEltsQuadratic(Elm)
	
  Nelt = size(Elm, 1);
	NgeoElm = zeros(Nelt*4, 3);
	cnt = 0;
  for ie = 1:Nelt

    nv = Elm(ie, :)';

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(1, 1);
		NgeoElm(cnt, 2) = nv(4, 1);
		NgeoElm(cnt, 3) = nv(6, 1);


		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(2, 1);
		NgeoElm(cnt, 2) = nv(5, 1);
		NgeoElm(cnt, 3) = nv(4, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(3, 1);
		NgeoElm(cnt, 2) = nv(6, 1);
		NgeoElm(cnt, 3) = nv(5, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(4, 1);
		NgeoElm(cnt, 2) = nv(5, 1);
		NgeoElm(cnt, 3) = nv(6, 1);

  end
		

return

function NgeoElm = subCompNewEltsCubic(Elm)

  Nelt = size(Elm, 1);
	NgeoElm = zeros(Nelt*9, 3);
	cnt = 0;
  for ie = 1:Nelt

		nv = Elm(ie, :)';

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(1, 1);
		NgeoElm(cnt, 2) = nv(4, 1);
		NgeoElm(cnt, 3) = nv(9, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(2, 1);
		NgeoElm(cnt, 2) = nv(6, 1);
		NgeoElm(cnt, 3) = nv(5, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(3, 1);
		NgeoElm(cnt, 2) = nv(8, 1);
		NgeoElm(cnt, 3) = nv(7, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(4, 1);
		NgeoElm(cnt, 2) = nv(5, 1);
		NgeoElm(cnt, 3) = nv(10, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(6, 1);
		NgeoElm(cnt, 2) = nv(7, 1);
		NgeoElm(cnt, 3) = nv(10, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(8, 1);
		NgeoElm(cnt, 2) = nv(9, 1);
		NgeoElm(cnt, 3) = nv(10, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(4, 1);
		NgeoElm(cnt, 2) = nv(10, 1);
		NgeoElm(cnt, 3) = nv(9, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(5, 1);
		NgeoElm(cnt, 2) = nv(6, 1);
		NgeoElm(cnt, 3) = nv(10, 1);

		cnt = cnt + 1;
		NgeoElm(cnt, 1) = nv(7, 1);
		NgeoElm(cnt, 2) = nv(8, 1);
		NgeoElm(cnt, 3) = nv(10, 1);

	end

return

