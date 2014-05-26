%%%%%%%%%%%%%%%%
% input values %
%%%%%%%%%%%%%%%%
% Ndm: geometry information
% NobN: Number of elements in the outer boundary
% NobD: Number of elements in the inner boundary
% BdyN: Outer boundary points before ordering (index)
% BdyD: Inner boundary points before ordering (index)
% SolBD: Solution values on outer boundary 
% BdyDP: Potential values on inner boundary
% Elm: connectivity information
%%%%%%%%%%%%%%%%%
% output values %
%%%%%%%%%%%%%%%%%
% nBdyN: reordered outer boundary values
% outerNodes: reordered outer boundary nodes

function [outerValues, innerValues] = ...
  fePostProcOrdering(Ndm, NobN, NobD, BdyN, BdyD, SolBd, BdyDP, Elm)


  %plot the geometry of the nodes on the exterior of the surface
  figure(234); hold on;
  %first plot the nodes of the whole surface in black
  for i = 1:length(Ndm)
    plot(Ndm(i,2),Ndm(i,3),'.k')
  end
  %then plot the outer boundary nodes
  for i = 1:NobN
    plot(Ndm(BdyN(i),2),Ndm(BdyN(i),3),'.r');
  end

  Nnd = size(Ndm, 1);
  NbdN = size(BdyN, 1);
  NbdD = size(BdyD, 1);
  Nelt = size(Elm, 1);
  Nnde = size(Elm, 2);	% 3:Linear, 6:Quadratic, 10:Cubic
  
  % Original boundary reallocation
  switch Nnde
   case 3  % Linear
    Odr = 1;
   case 6  % Quadratic
    Odr = 2;
   case 10 % Cubic
    Odr = 3;
  end % of switch

  
  % Givn nodes(geometric/interpolant) fill NdEltCell with elements
  % and order in each of them  
  NdEltCell = cell(Nnd, 2); % {Element No.}{th on the Element}
  % Map2New = zeros(NndF, 1);  % Map from old node to new extended node
  % Map2Old = zeros(NndF, 1);  % Map from new node to old extended node

  for ie = 1:Nelt
    for jt = 1:Nnde
      ndi = Elm(ie, jt);
      
      % Set Elt/th on each node
      NdEltCell{ndi, 1} = [NdEltCell{ndi, 1}; ie];
      NdEltCell{ndi, 2} = [NdEltCell{ndi, 2}; jt];
    end
  end

  % Given Element Set, find neighboring 
  [outerNodes, outerValues] = subPostProcOrdering(NobN, Nnde, Odr, ...
						  BdyN, SolBd, Elm, ...
						  NdEltCell);
  figure(444);
  subplot(2,1,1);
  plot(SolBd);
  subplot(2,1,2);
  plot(outerValues);
  
  %plot the geometry of the nodes on the exterior of the surface
  figure(235); hold on;
  %first plot the nodes of the whole surface in black
  for i = 1:length(Ndm)
    plot(Ndm(i,2),Ndm(i,3),'.k')
  end
  %then plot the outer boundary nodes
  for i = 1:length(outerValues)
    plot(Ndm(outerNodes(i),2),Ndm(outerNodes(i),3),'.r');
  end

  
  str = num2str(Odr);
  title(cat(2, 'Solution on outer boundary in approximation order', str));
  
       
       
  [innerNodes, innerValues] = subPostProcOrdering(NobD, Nnde, Odr, BdyD,...
  BdyDP, Elm, NdEltCell);
  
  figure(888);
  subplot(2,1,1);
  plot(BdyDP);
  subplot(2,1,2);
  plot(outerNodes);
  str = num2str(Odr);
  title(cat(2, 'Given potential with interpolation of order', str));
  
     
  
return
  

  


function [nBdy, nBdyP] = subPostProcOrdering(Nob, Nnde, Odr, Bdy, ...
					     Soln, Elm, NdEltCell) 

% gBdy: Geometric Bdy node indices
% gBuf: Buffer starts
gBdy = Bdy(1:Nob, 1);
mxV = max(gBdy);
mnV = min(gBdy);
offset = mnV - 1;
gBuf = zeros(mxV - offset, 1);


% Moving geometric nodes
% Filling gBuf with new poisition on gBuf(Bdy(i, 1), 1)
nGnode = zeros(Nob, 1);
Nbd = size(Bdy, 1);

nBdy = zeros(Nbd, 1);
nBdyP = zeros(Nbd, 1);

for i = 1:Nob
  in = Odr * i - (Odr - 1);  % Mapping from one-by-one 
			     % increaseing to order-by-order increasing
  obd = Bdy(i, 1);
  nBdy(in, 1) = obd;
  nBdyP(in, 1) = Soln(i, 1);
  gBuf(obd - offset, 1) = in;
end % i


nb = size(Bdy, 1);
for j = Nob+1:nb
  
  bdNd = Bdy(j);
  eT = NdEltCell{bdNd, 1};	% Element having node j
  tH = NdEltCell{bdNd, 2};	% tH th position in Element
  
  tH = tH - 3;                  % Since triangle has 3 geo-nodes
  tH = tH + (Odr - 2);          % To derive right preceeding node
  
  pnodeidx = floor(tH/(Odr-1)); % pnode is preceeding to bdNd.
  pos = rem(tH, Odr-1);         % pos == 1 (Quadratic), 1,2(Cubic)
  
  pnode = Elm(eT, pnodeidx);
  pnodepos = gBuf(pnode - offset, 1);
  
  newpos = pnodepos + pos + 1;
  nBdy(newpos, 1) = bdNd;
  nBdyP(newpos, 1) = Soln(j, 1);
  
end % j

return




%function NdvalOld = fePostProcOrdering(Nnd, NsolIn, NsolBd, BdyDP, Map2New)
%NdvalOld = zeros(Nnd, 1);
%AllNval = zeros(Nnd, 1);
%AllNval = [NsolIn;NsolBd;BdyDP];

%for k = 1:Nnd
%  NdvalOld(k, 1) = AllNval(Map2New(k, 2), 1);
%end

%return
