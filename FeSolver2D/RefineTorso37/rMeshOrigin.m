% Gen data 2D - new version for h-test

% 4 parts
% - Bottom bar
% - Left Axis
% - Right Axis
% - Roof

% diV: Number of bisections on a unit interval


function [NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF] = rMeshOrigin(dataset, drawFlag)

addpath 'mat';
addpath 'data';


file_nodes = cat(2, dataset, '.pts');
file_tries = cat(2, dataset, '.tris');

NbdN = 53;                  % Neumann Bdy   
NbdD = 52;                  % Dirichlet Bdy


NdmF = load(file_nodes);     % Nnd = 659 'utahTorso34';      % Number of Totoal Nodes 16/48/180
                            % Nnd = 689 'utahTorso37_g';


Nnd = size(NdmF, 1)

ElmF = load(file_tries);

% % Draw dataset
if drawFlag == 1
    figno = 5555;
    subPlotData(figno, Nnd, NdmF, ElmF);
end

[NdmF, ElmF] = subValidateConnectivity(NdmF, ElmF);

Nnd = size(NdmF, 1);

Nint = Nnd - NbdN - NbdD;   % Interior Nodes

NelF = size(ElmF, 1)

for i = 1:NelF
    n1 = ElmF(i, 1);
    n2 = ElmF(i, 2);
    n3 = ElmF(i, 3);
    
    ElmF(i, 1) = n1;
    ElmF(i, 2) = n3;
    ElmF(i, 3) = n2;
end

NdMat = sparse(Nnd, 2);     % [X, Y] Coordinate set of all nodes

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

InitBdDID = 638;

for i = 1:NbdD
    BdyD(i, 1) = InitBdDID + i - 1;
end
BdyD;

% Coordinate Data on Neumann boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BdyN = zeros(NbdN, 1);

InitBdNID = 1;
for i = 1:NbdN
    BdyN(i, 1) = InitBdNID + i - 1;
end


% Coordinate Data in Interior of Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nint = Nnd - NbdD - NbdN;

IntNd = zeros(Nint, 1);

cnt = 0;
for i = 1:Nnd
    if i > NbdN                    % i is not in Neumann Boundary
        if i < InitBdDID | i > InitBdDID + NbdD - 1 % i is not in Dirichlet Boundary
            cnt = cnt+1;
            IntNd(cnt, 1) = i;
        end
    end
end

IntNd;


% Begin:Testing to All D-bdy case
BdyD = [BdyN;BdyD];
BdyN = [];
NbdD = size(BdyD, 1);
NbdN = size(BdyN, 1);
% End:Testing to All D-bdy case


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
Ncd = Nnd;
NdcondF = zeros(Ncd, 4);
  for i = 1:Ncd
    cx = NdMat(i, 1);
    cy = NdMat(i, 2);
    
    % \sigma(x,y) =  %%%%%%%% Ex: Function 0 
    condmat =  feval('feUserConduct', cx, cy);
    NdcondF(i, 1) = condmat(1,1,1);
    NdcondF(i, 2) = condmat(1,1,2);
    NdcondF(i, 3) = condmat(1,1,3);
    NdcondF(i, 4) = condmat(1,1,4);
  end

% Save above data into files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ttmp = cputime;
 % Coordinates
 save mat/data_pts.mat NdMat % -ASCII
 
 % Elements
 save mat/data_elt.mat Elm % -ASCII
 
 % Dirichlet Boundary Coordinate Indices
 save mat/data_bdyD.mat BdyD  %-ASCII
 
 % Neumann Boundary Coordinate Indices
 save mat/data_bdyN.mat BdyN  %-ASCII
 
 % Interior Coordinate Indices
 save mat/data_intNd.mat IntNd  %-ASCII
 
 % Boundary Potential on Dirichlet Boundary
 save mat/data_bdyDP.mat BdyDP %-ASCII
 
 % Conductivity data for each node
 save mat/data_condc.mat NdcondF %-ASCII

return

function subPlotData(figno, Nnd, NdmF, ElmF) 

Xvec = NdmF(:,1)';
Yvec = NdmF(:,2)';

NelF = size(ElmF, 1);


  figure(figno);
    triplot(ElmF, Xvec, Yvec, 'Color', [0.2 0.5 0.5]); hold on;

%     Plot global nodes
    for i = 1:Nnd
        str = num2str(i);
        text(Xvec(1, i), Yvec(1, i), str, 'Color', [0.7 0.2 0], 'Fontsize', 13);
    end

    hold on;
    
%     Plot Element No.
    for i = 1:NelF
        one = ElmF(i,1); two = ElmF(i,2); thr = ElmF(i,3);
        xp = (Xvec(1, one) + Xvec(1, two) + Xvec(1, thr))/3;
        yp = (Yvec(1, one) + Yvec(1, two) + Yvec(1, thr))/3;
        str = num2str(i);
        text(xp, yp, str, 'Color', [0 0.5 0], 'Fontsize', 14);
    end

    hold on;
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
