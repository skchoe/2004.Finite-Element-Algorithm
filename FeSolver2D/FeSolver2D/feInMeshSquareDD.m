% Gen data 2D - new version for h-test

% 4 parts
% - Bottom bar
% - Left Axis
% - Right Axis
% - Roof

% diV: Number of bisections on a unit interval


function [NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF] = ...
         feInMeshSquareDD(leftE, rightE, diV)

  if diV == 1      Nint = 0; % Interior Nodes
  elseif diV == 2  Nint = 1;
  else             
      Nint = (diV-1)^2;
  end

  Nint = (diV-1)*(diV-1);
  NbdN = 0;   % Neumann Bdy   
  
  NbdD = 2*(diV + 1) + 2*(diV - 1);               % Dirichlet Bdy

  Nnd = Nint + NbdN + NbdD;               % Number of Totoal Nodes 16/48/180

  NdMat = zeros(Nnd, 2);                    % [X, Y] Coordinate set of all nodes


  % Fill in NdMat (Array of Nodes(Coordinates)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cnt = 0;
 
  Nx = diV + 1; % Number of horizontal/vertical nodes
  Ny = Nx;


  dxy = (rightE-leftE)/diV;

  for j = 1:Ny

    dy = leftE + (j-1) * dxy;
    for i = 1:Nx
        dx = leftE + (i-1) * dxy;
        cnt = cnt+1;
        % X,Y input
        NdMat(cnt, 1) = dx;
        NdMat(cnt, 2) = dy;
    end
    
  end

  NdMat;

  % Fill in Elm (Array of Elements) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Nelt = diV^2 * 2;     % Number of Elements
  Elm = zeros(Nelt, 3);   % [Node1, Node2, Node3]'

  cntE = 0;
  for j = 1:diV

    cntN1 = (j-1)*Nx;
    cntN2 = j*Nx;

    for i = 1:diV
        cntE = cntE + 1;
        Elm(cntE, 1) = cntN1 + i;
        Elm(cntE, 2) = cntN2 + i + 1;
        Elm(cntE, 3) = cntN2 + i;
            
        cntE = cntE + 1;
        Elm(cntE, 1) = cntN1 + i;
        Elm(cntE, 2) = cntN1 + i + 1;
        Elm(cntE, 3) = cntN2 + i + 1;
    end
  end

Elm;

% Coordinate Data on Neumann boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BdyN = [];

% Coordinate Data on Dirichlet boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BdyD = zeros(NbdD, 1);

startN = 0;
cnt = 0;
for i = 1:diV + 1
    cnt = cnt + 1;
    BdyD(cnt, 1) = startN + i;
end

startN = startN + diV + 1;

for j = 1:diV - 1
    cnt = cnt + 1;
    BdyD(cnt, 1) = startN + 1;

    startN = startN + diV + 1;
    cnt = cnt + 1;
    BdyD(cnt, 1) =  startN;
    
end

for i = 1:diV + 1
    cnt = cnt + 1;
    BdyD(cnt, 1) = startN + i;
end
BdyD;


% Coordinate Data in Interior of Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nint = Nnd - NbdD - NbdN;

IntNd = zeros(Nint, 1);
intidx = 0;
startI = diV + 2;

if Nint ~= 0
  for iy = 1:diV-1
    for ix = 1:diV-1
        
      intidx = intidx + 1;
      IntNd(intidx, 1) = startI + ix;
        
    end
    startI = startI + diV + 1;
  end
end
  
IntNd;
    
  
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

    BdyDP(i, 1) = feInSolution(cx, cy);
  end
BdyDP;





% Conductivity on each node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output:   Ndcond = Array of Given Conductivity on each element
% input :   Ndcond = Same of output 
%           Nbp = Size of Dirichlet Boundary
%           BdyD = Array of Indices of Dirichlet boundary in Coordinate Set

% Conductivity may be given by  \sigma(ei, x,y) = xy
% or                            \sigma(ei, x,y) = 1/(xy) if this is well defined, x, y is in element ei
Ncd = Nelt;
NdcondF = zeros(3, 1, Ncd, 4);     % NdcondF contains
  for ei = 1:Ncd
    NdinE = Elm(ei,:);          % [1, 5, 4]
    cxv = NdMat(NdinE, 1);      % [-1 0 -1]'
    cyv = NdMat(NdinE, 2);      % [-1 0  0]'

    % \sigma(x,y) =  %%%%%%%% Ex: Function 0;
    condmat =  feInConduct(ei, cxv, cyv);
    for j = 1:3
        NdcondF(j, 1, ei, 1) = condmat(j,1,1);
        NdcondF(j, 1, ei, 2) = condmat(j,1,2);
        NdcondF(j, 1, ei, 3) = condmat(j,1,3);
        NdcondF(j, 1, ei, 4) = condmat(j,1,4);
    end

  end
  
NdcondF;



% Save above data into files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Coordinates
%save mat/data_pts.mat NdMat %-ASCII
%
%% Elements
%save mat/data_elt.mat Elm  %-ASCII
%
%% Dirichlet Boundary Coordinate Indices
%save mat/data_bdyD.mat BdyD  %-ASCII
%
%% Neumann Boundary Coordinate Indices
%save mat/data_bdyN.mat BdyN  %-ASCII
%
%% Interior Coordinate Indices
%save mat/data_intNd.mat IntNd  %-ASCII
%
%% Boundary Potential on Dirichlet Boundary
%save mat/data_bdyDP.mat BdyDP %-ASCII
%
%% Conductivity data for each node
%save mat/data_condc.mat NdcondF %-ASCII
%

return
