% Gen data 2D - new version for h-test
% utahTorso34 Info
% Total geometric node = 659 (But 659th node is never referenced)
% Nodes which is unreferenced = 1;
% Number of Elements = 1203
% Number of Neumann Boundary : (53) 1-51, 580, 650
% Number of Dirichlet Boundary : (60) 198 - 257
% Number of Interior nodes : (545) {52, ...,197,258,...579,581,...,649,651,...658} 


% 4 parts
% - Bottom bar
% - Left Axis
% - Right Axis
% - Roof

% diV: Number of bisections on a unit interval


function [NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF] = rMeshOrigin34(dataset, drawFlag)

addpath 'mat';
addpath 'data';
addpath '../FeSolver2D';


file_nodes = cat(2, dataset, '.pts');
file_tries = cat(2, dataset, '.tris');

NdmF = load(file_nodes);      %  Nnd = 659;      % Number of Totoal Nodes 16/48/180

Nnd = size(NdmF, 1);

ElmF = load(file_tries);

% % Draw dataset
if drawFlag == 1
    figno = 55;
    subPlotData(figno, Nnd, NdmF, ElmF);
end

[NdmF, ElmF] = subValidateConnectivity(NdmF, ElmF);

Nnd = size(NdmF, 1);
NelF = size(ElmF, 1);

for i = 1:NelF
    n1 = ElmF(i, 1);
    n2 = ElmF(i, 2);
    n3 = ElmF(i, 3);
    
    ElmF(i, 1) = n1;
    ElmF(i, 2) = n3;
    ElmF(i, 3) = n2;
end

NbdN = 53;                  % Neumann Bdy   
NbdD = 60;                  % Dirichlet Bdy
Nint = Nnd - NbdN - NbdD;   % Interior Nodes

NdMat = sparse(Nnd, 2);                    % [X, Y] Coordinate set of all nodes

% % % GEN_DATA_GENSTART = 0
Ttmp = cputime;
% Fill in NdMat (Array of Nodes(Coordinates)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NdMat = NdmF;

% % % GEN_NODESET_GEND = cputime - Ttmp

% Fill in NdMat (Array of Nodes(Coordinates)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nelt = NelF;     % Number of Elements
Elm = ElmF;   % [Node1, Node2, Node3]'%'

% % % % GEN_ELEMENTSET_GEND = cputime - Ttmp
% % % % Ttmp = cputime;

% Coordinate Data on Dirichlet boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BdyD = zeros(NbdD, 1);

InitBdDID = 198;

for i = 1:NbdD
    BdyD(i, 1) = InitBdDID + i - 1;
end


% Coordinate Data on Neumann boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BdyN = zeros(NbdN, 1);

InitBdNID = 1;
for i = 1:NbdN-2
    BdyN(i, 1) = InitBdNID + i - 1;
end

% Exception
spIdx1 = 580;
spIdx2 = 650;
BdyN(NbdN-1, 1) = spIdx1;
BdyN(NbdN, 1) = spIdx2;


% Coordinate Data in Interior of Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nint = Nnd - NbdD - NbdN;

IntNd = zeros(Nint, 1);

cnt = 0;

for i = 1:Nnd
    if i > NbdN-2 & i~= spIdx1 & i~= spIdx2   % i is not in Neumann Boundary
        if i < InitBdDID | i > InitBdDID + NbdD - 1 % i is not in Dirichlet Boundary
            cnt = cnt+1;
            IntNd(cnt, 1) = i;
        end
    end
end
IntNd;


% Testing to All D-bdy case
% BdyD = [BdyN;BdyD];
% BdyN = [];
NbdD = size(BdyD, 1);
NbdN = size(BdyN, 1);

% Interior = IntNd'
% BdyDirichlet = BdyD'
% BdyNeumann = BdyN'
  
% % % % GEN_BDYINDICES_SET_GEND = cputime - Ttmp


% Potential on Dirichlet boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output:   BdyDP = Array of Given Potential on Dirichlet Boundary
% input :   BdyDP = Same of output 
%           Nbp = Size of Dirichlet Boundary
%           BdyD = Array of Indices of Dirichlet boundary in Coordinate Set

Nbp = NbdD;
BdyDP = zeros(Nbp, 1);

  for i = 1:Nbp
    idx = BdyD(i, 1);
    cx = NdMat(idx, 1);
    cy = NdMat(idx, 2);

    BdyDP(i, 1) = feval('feUserSolution', cx, cy);
  end
BdyDP;

% Conductivity on each node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output:   Ndcond = Array of Given Conductivity on each node
% input :   Ndcond = Same of output 
%           Nbp = Size of Dirichlet Boundary
%           BdyD = Array of Indices of Dirichlet boundary in Coordinate Set

% Conductivity may be given by  \sigma(x,y) = xy
% or                            \sigma(x,y) = 1/(xy) if this is well defined

  Ncd = Nelt;
  NdcondF = zeros(3, 1, Ncd, 4);     % NdcondF contains
  for ei = 1:Ncd
    NdinE = Elm(ei,:);          % [1, 5, 4]
    cxv = NdMat(NdinE, 1);      % [-1 0 -1]'
    cyv = NdMat(NdinE, 2);      % [-1 0  0]'

    % \sigma(x,y) =  %%%%%%%% Ex: Function 0;
    condmat =  feval('feUserConduct', ei, cxv, cyv);
    for j = 1:3
        NdcondF(j, 1, ei, 1) = condmat(j,1,1);
        NdcondF(j, 1, ei, 2) = condmat(j,1,2);
        NdcondF(j, 1, ei, 3) = condmat(j,1,3);
        NdcondF(j, 1, ei, 4) = condmat(j,1,4);
    end
    
  end


% % % % GEM_FILES_WRITTEN = cputime - Ttmp

return

function subPlotData(figno, Nnd, NdmF, ElmF) 

Xvec = NdmF(:,1)';
Yvec = NdmF(:,2)';

NelF = size(ElmF, 1);


  figure(figno);
    triplot(ElmF, Xvec, Yvec, 'Color', [0.9 0.9 0.9]); hold on;

%     Plot global nodes
    for i = 1:Nnd
        str = num2str(i);
        text(Xvec(1, i), Yvec(1, i), str, 'Color', [0.5 0 0], 'Fontsize', 12);
    end

    hold on;
%     
% %     Plot Element No.
%     for i = 1:NelF
%         one = ElmF(i,1); two = ElmF(i,2); thr = ElmF(i,3);
%         xp = (Xvec(1, one) + Xvec(1, two) + Xvec(1, thr))/3;
%         yp = (Yvec(1, one) + Yvec(1, two) + Yvec(1, thr))/3;
%         str = num2str(i);
%         text(xp, yp, str, 'Color', [0 0.5 0], 'Fontsize', 7);
%     end
% 
%     hold on;
return

function [newNdmF, newElmF] = subValidateConnectivity(NdmF, ElmF)

  Nnd = size(NdmF, 1);
  Nelt = size(ElmF, 1);

  newElmF = zeros(Nelt, 3);
  
  
  ndTgSet = zeros(Nnd, 1);

  for ie = 1:Nelt
    for in = 1:3
        node = ElmF(ie, in);
        if ndTgSet(node, 1) == 0
            ndTgSet(node, 1) = 1;
        end
    end
  end
                
  newNdmF = [];
  Map2New = [];

  decrf = 0;
  newid = 0;
  for in = 1:Nnd
      newid = newid + 1;
      
      if ndTgSet(in, 1) == 0
          decrf = decrf + 1;
          Map2New(in, 1) = -1;  % This means this node should be excluded.
      else
          newNdmF = [newNdmF; NdmF(in, :)];
          Map2New(in, 1) = in - decrf;
      end
  end
  
  
  for ie = 1:Nelt
      for j = 1:3
          eid = ElmF(ie, j);
          newElmF(ie, j) = Map2New(eid, 1);
      end
  end
    
return
