% Gen data 2D - new version for h-test

% 4 parts
% - Bottom bar
% - Left Axis
% - Right Axis
% - Roof

% diV: Number of bisections on a unit interval


function [NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF] = ...
    feInMeshSqulus(diV)

addpath 'mat'
divisions = 2^(diV);
divisions2 = 2*divisions;    
divisions3 = 3*divisions;

Nint = (divisions - 1)^2 * 8 + (divisions - 1) * 8; % Interior Nodes
NbdN = 2*(divisions3 + 1) + 2*(divisions3 - 1);     % Neumann Bdy   
NbdD = (divisions - 1) * 4 + 4;                     % Dirichlet Bdy

Nnd = Nint + NbdN + NbdD;                 % Number of Totoal Nodes 16/48/180


NdMat = sparse(Nnd, 2);                    % [X, Y] Coordinate set of all nodes

% % % GEN_DATA_GENSTART = 0
Ttmp = cputime;
% Fill in NdMat (Array of Nodes(Coordinates)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0;

Nx = divisions3 + 1; % Number of horizontal/vertical nodes
Ny = Nx;

dxy = (1/2)^diV;
for j = 1:Ny

    dy = (j-1) * dxy;
    if j <= divisions + 1         % Lower part 
        for i = 1:Nx
            dx = (i-1) * dxy;
            cnt = cnt+1;
            % X,Y input
            NdMat(cnt, 1) = dx;
            NdMat(cnt, 2) = dy;
        end
        
    elseif j < divisions2 + 1   % Middle part
        for i = 1:divisions + 1
            dx = (i-1) * dxy;
            cnt = cnt+1;
            NdMat(cnt, 1) = dx;
            NdMat(cnt, 2) = dy;
        end
        for i = divisions2 + 1: divisions3 + 1
            dx = (i-1) * dxy;
            cnt = cnt + 1;
            NdMat(cnt, 1) = dx;
            NdMat(cnt, 2) = dy;
        end
    else                    % Upper part
        for i = 1:Nx
            dx = (i-1) * dxy;
            cnt = cnt+1;
            % X,Y input
            NdMat(cnt, 1) = dx;
            NdMat(cnt, 2) = dy;
        end       
    end
end

% % % GEN_NODESET_GEND = cputime - Ttmp



% Fill in NdMat (Array of Nodes(Coordinates)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nelt = 2^(2*diV) * 2 * 8;     % Number of Elements
Elm = zeros(Nelt, 3);   % [Node1, Node2, Node3]'%'

cntE = 0;
cntN1 = 0;
cntN2 = 0;
for j = 1: 3*(divisions)

    if j <= divisions         % Lower part 
        cntN1 = (j-1)*Nx;
        cntN2 = j*Nx;
        
        for i = 1:divisions3
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN2 + i + 1;
            Elm(cntE, 3) = cntN2 + i;
            
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN1 + i + 1;
            Elm(cntE, 3) = cntN2 + i + 1;
        end
        
    elseif j == divisions + 1

        cntN1 = divisions * Nx;
        cntN2 = cntN1 + Nx;

        for i = 1:divisions
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN2 + i + 1;
            Elm(cntE, 3) = cntN2 + i;

            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN1 + i + 1;
            Elm(cntE, 3) = cntN2 + i + 1;
        end

        cntN1 = cntN1 + divisions2;
        cntN2 = cntN2 + divisions + 1;
        for i = 1:divisions
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN2 + i + 1;
            Elm(cntE, 3) = cntN2 + i;
            
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN1 + i + 1;
            Elm(cntE, 3) = cntN2 + i + 1;
        end

    elseif j < divisions2

        cntN1 = cntN1 + divisions + 1;
        cntN2 = cntN2 + divisions + 1;

        for i = 1:divisions
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN2 + i + 1;
            Elm(cntE, 3) = cntN2 + i;

            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN1 + i + 1;
            Elm(cntE, 3) = cntN2 + i + 1;
        end

        cntN1 = cntN1 + divisions + 1;
        cntN2 = cntN2 + divisions + 1;
        for i = 1:divisions
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN2 + i + 1;
            Elm(cntE, 3) = cntN2 + i;
            
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN1 + i + 1;
            Elm(cntE, 3) = cntN2 + i + 1;
        end
        
    elseif j == divisions2

        cntN1 = cntN1 + divisions + 1;
        cntN2 = cntN2 + divisions + 1;

        for i = 1:divisions
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN2 + i + 1;
            Elm(cntE, 3) = cntN2 + i;

            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN1 + i + 1;
            Elm(cntE, 3) = cntN2 + i + 1;
        end

        cntN1 = cntN1 + divisions + 1;
        cntN2 = cntN2 + divisions2;
        for i = 1:divisions
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN2 + i + 1;
            Elm(cntE, 3) = cntN2 + i;
            
            cntE = cntE + 1;
            Elm(cntE, 1) = cntN1 + i;
            Elm(cntE, 2) = cntN1 + i + 1;
            Elm(cntE, 3) = cntN2 + i + 1;
        end
        
    else
        
        if j==divisions2+1
            cntN1 = cntN1 + divisions + 1;
            cntN2 = cntN2 + divisions + 1;
        else
            cntN1 = cntN1 + divisions3 + 1;
            cntN2 = cntN2 + divisions3 + 1;
        end

        for i = 1:divisions3
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
end

Elm;


% % % % GEN_ELEMENTSET_GEND = cputime - Ttmp
% % % % Ttmp = cputime;

% Coordinate Data on Dirichlet boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BdyD = zeros(NbdD, 1);

startN = (divisions3 + 1)*divisions + divisions;
cnt = 0;
for i = 1:divisions + 1
    cnt = cnt + 1;
    BdyD(cnt, 1) = startN + i;
end

cntJ = divisions - 1;
startN = startN + divisions3+1;

for j = 1:cntJ
    cnt = cnt + 1;
    BdyD(cnt, 1) = startN + 1;
    cnt = cnt + 1;
    BdyD(cnt, 1) = BdyD(cnt-1, 1) + 1;
    startN = startN + divisions2 + 2;
end

for i = 1:divisions + 1
    cnt = cnt + 1;
    BdyD(cnt, 1) = startN + i;
end
BdyD;


% Coordinate Data on Neumann boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BdyN = zeros(NbdN, 1);
cnt = 0;
for i = 1:divisions3 + 1
    cnt = cnt + 1;
    BdyN(cnt, 1) = i;
end

startN = i;

for j = 1:divisions
    cnt = cnt+1;
    BdyN(cnt, 1) = startN + 1;
    cnt = cnt+1;
    BdyN(cnt, 1) = BdyN(cnt-1, 1) + divisions3;
    startN = BdyN(cnt, 1);
end

for k = 1:divisions - 1
    cnt = cnt + 1;
    BdyN(cnt, 1) = startN + 1;
    cnt = cnt + 1;
    BdyN(cnt, 1) = BdyN(cnt-1, 1) + divisions2 + 1;
    startN = BdyN(cnt, 1);
end

for l = 1:divisions
    cnt = cnt + 1;
    BdyN(cnt, 1) = startN + 1;
    cnt = cnt+1;
    BdyN(cnt, 1) = BdyN(cnt-1, 1) + divisions3;
    startN = BdyN(cnt, 1);
end

 for m = 1:divisions3 + 1
    cnt = cnt + 1;
    BdyN(cnt, 1) = startN + m;
end



% Coordinate Data in Interior of Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nint = Nnd - NbdD - NbdN;

IntNd = zeros(Nint, 1);

if Nint ~= 0

  flg = 0;
  intidx = 0;

  for i = 1:Nnd
    if flg == 0
        for k=1:NbdN
            if i == BdyN(k, 1)
              flg = 1;
              break;
            end
        end
    end
    
    if flg == 0
        for j=1:NbdD
            if i == BdyD(j, 1)
                flg = 1;
              break;
            end
        end
    end

    if flg == 0
      intidx = intidx + 1;
      IntNd(intidx, 1) = i;
    end
    
    flg = 0;
    
  end

end

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

    BdyDP(i, 1) = feInSolution(cx, cy);
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
NdcondF = zeros(3, 1, Ncd, 4);
  for i = 1:Ncd
		NdinE = Elm(ei,:);
    cxv = NdMat(NdinE, 1);
    cyv = NdMat(NdinE, 2);
    
    % \sigma(x,y) =  %%%%%%%% Ex: Function 0
    condmat =  feInConduct(ei, cxv, cyv);
		for j = 1:3
      NdcondF(j, 1, ei, 1) = condmat(j,1,1);
      NdcondF(j, 1, ei, 2) = condmat(j,1,2);
      NdcondF(j, 1, ei, 3) = condmat(j,1,3);
      NdcondF(j, 1, ei, 4) = condmat(j,1,4);
    end

 end

% Save above data into files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ttmp = cputime;
% Coordinates
%save mat/data_pts.mat NdMat % -ASCII

%% Elements
%save mat/data_elt.mat Elm % -ASCII

% Dirichlet Boundary Coordinate Indices
%save mat/data_bdyD.mat BdyD  %-ASCII

% Neumann Boundary Coordinate Indices
%save mat/data_bdyN.mat BdyN  %-ASCII

%% Interior Coordinate Indices
%save mat/data_intNd.mat IntNd  %-ASCII

%% Boundary Potential on Dirichlet Boundary
%save mat/data_bdyDP.mat BdyDP %-ASCII

%% Conductivity data for each node
%save mat/data_condc.mat NdcondF %-ASCII

% % % % GEM_FILES_WRITTEN = cputime - Ttmp

return
