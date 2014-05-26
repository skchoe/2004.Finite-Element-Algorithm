% Ndm, Elm, BdyD, BdyN, IntNd, BdyDP, Elcond 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ndm, Elm, BdyN, BdyD, IntNd, BdyDP, Elcond, Ndmap2New] = ...
    feDataProc(bsOdr, NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, Elcond, drFlag)
%------------------------------------------------------

  % NLcomp : Number of Nodes in Local Element
  switch bsOdr
  case 'L'
    NLcomp = 3;
  case 'Q'
    NLcomp = 6;
  case 'C'
    NLcomp = 10;
  end
    
% %     % Data loading and store in primitive structure
% % %     PPROC_PTS_ELT_LOADED = 0
    
  % Ndm:Nodes Coordinates initialization------------------
  NndF = size(NdMat, 1);
  NdmF = zeros(NndF, 4);

  for i = 1:NndF
    NdmF(i, 1) = i;
  end

  NdmF(:,2) = NdMat(:,1);
  NdmF(:,3) = NdMat(:,2);

  % PPROC_PTS_ELT_ALLOCATED = 0

  %------------------------------------------------------
  % ElmF:Elements initialization
  NeltF = size(Elm, 1);
  ElmF = zeros(NeltF, 3); % Elt size of geometric element = 3
  ElmF = Elm;
  
    
  % Neumann Boundary Indices setting
  BdyNF = BdyN;

  % Dirichlet Boundary Indices setting
  BdyDF = BdyD;

% % Interior Node indices setting
  IntNdF = IntNd;
  if size(IntNdF) == 0
    IntNdF = [];
  end

  % Dirichlet Boundary setting
  BdyDPF = BdyDP;

  % Conductivity on each nodes
  ElcondF = Elcond; % Variable name are same as stored object name

    
  % Generate node, ... by the value of interpolation order
  [xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xElcond] ...
    = subSetHigherNodes(bsOdr, NdmF, ElmF, IntNdF, BdyNF, BdyDF, BdyDPF, ElcondF);

  Ndmap2New = zeros(1,1);
   
  % Node, Element setting - Ordering routine for data 
  [Ndm, Elm, IntNd, BdyN, BdyD, BdyDP, Elcond, Ndmap2New] ...
    = subPreProcOrdering(NLcomp, xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xElcond);

  if drFlag == 1
    subPlotAppMesh(133, NLcomp, Ndm, Elm);
  end

return

% This function gets input from file data and generates new nodes depending on approximation error
function [xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xElcond] = ...
    subSetHigherNodes(bsOdr, NdmF, ElmF, IntNdF, BdyNF, BdyDF, BdyDPF, ElcondF)

  NndF = size(NdmF, 1);            % [ID, X, Y, U, NAdjEdges]'
  NelmF = size(ElmF, 1);           % [Node1, Node2, Node3]'

  Nint = size(IntNdF, 1);
  NbdN = size(BdyNF, 1);
  NbdD = size(BdyDF, 1);

  NbdDP = size(BdyDPF, 1);
  Ncond = size(ElcondF, 1);
% % % % % % % %   
% % % % % % % %   % Test code I'm working on!
% % % % % % % %   ElconElF = ssubSetAdjElement(NdmF, ElmF);
% % % % % % % % 
  % IntNdF, BdyNF, BdyDF, BdyDPF, NdcondF values are newly
  % generated from old values. 
  
  
  % NLcomp : Number of Nodes in Local Element
  switch bsOdr
    case 'L'
      [xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xElcond]...
          = ssubSetHnodesLinear(NdmF, ElmF, IntNdF, BdyNF, BdyDF, ...
				BdyDPF, ElcondF);   
    case 'Q'
      [xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xElcond]...
          = feSetHnodesQuadratic(NdmF, ElmF, IntNdF, BdyNF, BdyDF, ...
				 BdyDPF, ElcondF);     
    case 'C'
      [xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xElcond]...
          = feSetHnodesCubic(NdmF, ElmF, IntNdF, BdyNF, BdyDF, BdyDPF, ElcondF);
  end
  

return


function [xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xElcond]...
    = ssubSetHnodesLinear(NdmF, ElmF, IntNdF, BdyNF, BdyDF, BdyDPF, ElcondF)

  xNdm = NdmF;
  xElm = ElmF;
  xIntNd = IntNdF;
  xBdyN = BdyNF;
  xBdyD = BdyDF;
  xBdyDP = BdyDPF;
  xElcond = ElcondF;

return

function subPlotGeoMesh(figno, Ndm, Elm)

  Xvec = Ndm(:,1)';
  Yvec = Ndm(:,2)';

  % Global Node, Elt sturcture plotting
  figure(figno);
    triplot(Elm, Xvec, Yvec);   hold on;
    plotgnode(Xvec, Yvec);      hold on;
    ploteltno(Elm, Xvec, Yvec); hold on;

return

function subPlotAppMesh(figno, NLcomp, Ndm, Elm)

  Xvec = Ndm(:,2)';
  Yvec = Ndm(:,3)';
    
  Nelm = size(Elm, 1);
  vec1 = zeros(1,3);
  vec2 = zeros(1,3);


% Normal test in each element
%   for i = 1:Nelm
%       v1 = Elm(i, 1);
%       v2 = Elm(i, 2);
%       v3 = Elm(i, 3);
%       
%       vc11 = Ndm(v2, 2) - Ndm(v1, 2);
%       vc12 = Ndm(v2, 3) - Ndm(v1, 3);
%       vc21 = Ndm(v3, 2) - Ndm(v1, 2);
%       vc22 = Ndm(v3, 3) - Ndm(v1, 3);
%       vec1(1,1) = vc11; vec1(1,2) = vc12;
%       vec2(1,1) = vc21; vec2(1,2) = vc22;
%
%       crvec = cross(vec1, vec2)
%   end
%
%   Nelm = size(Elm, 1)

  figure(figno);
    triplot(Elm(:,1:3), Xvec, Yvec); hold on;
    plotgnode(Xvec, Yvec);      hold on;
    ploteltno(Elm, Xvec, Yvec); hold on;
    
return

%---------------------------------------------------------
% This function gets node info which is determined by approximation
% order and reorder it which we can distinguish it to interior/bdyN/bdyD.
function [newNdm, newElm, newIntNd, newBdyN, newBdyD, newBdyDP, ...
	  newElcond, Ndmap2New] = subPreProcOrdering(NLcomp, xNdm, ...
						  xElm, xIntNd, ...
						  xBdyN, xBdyD, xBdyDP, xElcond)

  Nnd = size(xNdm, 1);            % [ID, X, Y, U, NAdjEdges]'
  Nelm = size(xElm, 1);           % [Node1, Node2, Node3]'

  Nint = size(xIntNd, 1);
  NbdN = size(xBdyN, 1);
  NbdD = size(xBdyD, 1);

  NbdDP = size(xBdyDP, 1);
  
  newNdm = zeros(Nnd, 4);
  newElm = zeros(Nelm, NLcomp);

  %-------------------------------
  % Rebuilding Ndm/IntNd/BdyN/BdyD 
  % Setup the mapping table
  Ndmap2New = zeros(Nnd, 2);
  Ndmap2Old = zeros(Nnd, 2);

  ndorder = 0;

  if Nint > 0
    for i = 1:Nint

        ndorder = ndorder + 1;

% %         if xIntNd(i, 1) > SizeofnewNodeMat | ndorder > SizeofnewNodeMat
% %             ndorder
% %             xIntNd = xIntNd(i, 1)
% %         end
        newNdm(ndorder,:) = xNdm(xIntNd(i, 1), :);
        newNdm(ndorder, 1) = ndorder;
        newIntNd(i, 1) = ndorder;
        
        Ndmap2Old(ndorder, 1) = ndorder;
        Ndmap2Old(ndorder, 2) = xIntNd(i, 1);
        Ndmap2New(xIntNd(i, 1), 1) = xIntNd(i, 1);        
        Ndmap2New(xIntNd(i, 1), 2) = ndorder;

    end
  else
      newIntNd = [];
  end
  
  if NbdN ~= 0
      
    for j = 1:NbdN
      ndorder = ndorder + 1;
      
      newNdm(ndorder, :) = xNdm(xBdyN(j, 1), :);
      newNdm(ndorder, 1) = ndorder;
      newBdyN(j, 1) = ndorder;
      
      Ndmap2Old(ndorder, 1) = ndorder;
      Ndmap2Old(ndorder, 2) = xBdyN(j, 1);
      Ndmap2New(xBdyN(j, 1), 1) = xBdyN(j, 1);        
      Ndmap2New(xBdyN(j, 1), 2) = ndorder;
    end
  else
      newBdyN = [];
  end
  
  if NbdD ~= 0
    for k = 1:NbdD
      ndorder = ndorder + 1;
      
      newNdm(ndorder, :) = xNdm(xBdyD(k, 1), :);
      newNdm(ndorder, 1) = ndorder;
      newBdyD(k, 1) = ndorder;
      
      Ndmap2Old(ndorder, 1) = ndorder;
      Ndmap2Old(ndorder, 2) = xBdyD(k, 1);
      Ndmap2New(xBdyD(k, 1), 1) = xBdyD(k, 1);        
      Ndmap2New(xBdyD(k, 1), 2) = ndorder;
    end
    
  else    
    newBdyD = [];
  end

  %-------------------------------------------------------
  % Node indices update
  
  for l = 1:Nelm
    for m = 1:NLcomp
        newElm(l, m) = Ndmap2New(xElm(l, m), 2);
    end
  end
  
  %-------------------------------------------------------
  % Conductivity/Boundary Values order update 
  % - No change (Ordering change affects index, not value on the point)

  newElcond = xElcond;
  newBdyDP = xBdyDP;

return

%------------------------------------------------------
function ElconEl = ssubSetAdjElement(Ndm, Elm)

  Nnd = size(Ndm, 1);
  Nelt = size(Elm, 1);


  % Fill in edge list to each node where node is included.
  NShareNode = 0;

  for i = 1:Nelt
    for j = 1:3
        Nde = Elm(i, j);
        Ndm(Nde, 4) = Ndm(Nde, 4) + 1;
    end    
  end

  MxConNd = max(Ndm(:, 4));

  NdconEl = zeros(Nnd, MxConNd+1);      % [NumConn, ConElt1, ConElt2, ...]

  for i = 1:Nelt
    for j = 1:3                         % For each node of given element
        Nde = Elm(i, j);
        NdconEl(Nde, 1) = NdconEl(Nde, 1) + 1;
        offset = 1 + NdconEl(Nde, 1);   % Position where ConElti is stored
        NdconEl(Nde, offset) = i;       % Put this element to
                                        % NdconEl list of given
                                        % node 'Nde' 
    end    
  end

  ElconEl = zeros(Nelt, 3);             % [ConnNode1(eg1),
                                        % ConnNode2(eg2),
                                        % ConnNode3(eg3)] 

  %------------------------------------------------------
  % Check connectivity across boundary of elements
  conEl = 0;
  th = 0;
  % Nelt = 2;
  for i = 1:Nelt                        % For each element,
    for k = 1:3                         % For each edge of given element

        eg1 = k;
        if k==3 eg2 = 1; else eg2 = k+1; end

        Nd1 = Elm(i, eg1);
        Nd2 = Elm(i, eg2);

        Nel = NdconEl(Nd1, 1);          % Number of element containing Node1
        % Find Element having Node2 except elt i 
	[conEl, th] = sssubFindNodeIn(Nd2, i, Nel, NdconEl(Nd1,:), Elm); 
        ElconEl(i, k) = conEl;          % Put the element to
                                        % connecting neighbor of
                                        % elt i 
        
    end
  end

return

% output: conEl = output of ?th node in element discovered
% input : Nd = Node 
% eltin : original element
% Num   : number of valid connecting element
% Ndvec : vector of elements containing calling node
% Elm   : Table of Nodes for each element
function [conEl, th] = sssubFindNodeIn(Nd2, eltin, Num, Eltvec, Elm)

    conEl = 0;
    th = 0;
    
    for i = 1:Num                       % For each element

        elt = Eltvec(1, i+1);
        if elt ~= eltin
            for j = 1:3                 % Find elt has node Nd2
                node = Elm(elt, j);
                if Elm(elt, j) == Nd2
                    conEl = elt;
                    th = j;
                    break;
                end
            end
        end
        
    end
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Node numbers plotting
function plotgnode(Xvec, Yvec)

    Nnd = size(Xvec, 2);
    for i = 1:Nnd
        str = num2str(i);
        text(Xvec(1, i), Yvec(1, i), str, 'Color', [0.5 0 0], 'Fontsize', 8);
    end
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Element numbers plotting
function ploteltno(Elm, Xvec, Yvec)

    Nel = size(Elm, 1);
    for i = 1:Nel
        one = Elm(i,1); two = Elm(i,2); thr = Elm(i,3);
        xp = (Xvec(1, one) + Xvec(1, two) + Xvec(1, thr))/3;
        yp = (Yvec(1, one) + Yvec(1, two) + Yvec(1, thr))/3;
        str = num2str(i);
        text(xp, yp, str, 'Color', [0 0.5 0], 'Fontsize', 7);
    end

return

