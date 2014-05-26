% Gen data 2D - new version for h-test
% layer70_4 Info
% Total geometric node = 3859 
% Nodes which is unreferenced = 0
% Number of Elements = 6893
% Number of Neumann Boundary : (105) 1-105
% Number of Dirichlet Boundary : (720) (given in a file)
% Number of Interior nodes : 3034


% 4 parts
% - Bottom bar
% - Left Axis
% - Right Axis
% - Roof

% diV: Number of bisections on a unit interval


function [geometry, connectivities, outerHeart, outerTorso, interior, ...
	  heartPotential, conductivity] =  feMeshTorsoOrg340(dataset, ...
						  Elcond, drawFlag, ...
						  isFakeSol) 

addpath 'data';

plotMesh = 0;

file_nodes = cat(2, dataset, '.pts');
file_connectivity = cat(2, dataset, '.fac');
heartBoundaryNodes= load('boundary.pts');

geometry = load(file_nodes);      % numPoints = 3859;     
				  % Number of Total Nodes 
				  
numPoints = size(geometry, 1);

connectivities = load(file_connectivity);

[geometry, connectivities] = subValidateConnectivity(geometry, ...
						  connectivities); 

numPoints = size(geometry, 1);
numElements = size(connectivities, 1);

% reverse the order of the triangle connectivity
% from clockwise to counter-clockwise
connectivities = permute(connectivities, [1;3;2]);

numNeumann = 105;                                      % Neumann Bdy   
numDirichlet = 720;                                    % Dirichlet Bdy
numInterior = numPoints - numNeumann - numDirichlet;   % Interior Nodes

% Coordinate Data on Dirichlet boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outerHeart = heartBoundaryNodes;

% Coordinate Data on Neumann boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outerTorso = (1:numNeumann)';

% Coordinate Data in Interior of Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interior = zeros(numInterior, 1);

cnt = 1;
for i = 1:numPoints
  % i is not in the Neumann boundary
  if i > numNeumann
    % i is not in the Dirichlet boundary
    if isempty(find(outerHeart == i))
      interior(cnt, 1) = i; cnt = cnt+1;
    end
  end
end


if isFakeSol == 1
  outerHeart = [outerTorso;outerHeart];
  outerTorso = []; 
end

numDirichlet = size(outerHeart, 1); 
numNeumann = size(outerTorso, 1);


% Potential on Dirichlet boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output:   heartPotential = Array of Given Potential on Dirichlet Boundary
% input :   heartPotential = Same of output 
%           outerHeart = Array of Indices of Dirichlet boundary in
%           Coordinate Set 

heartPotential = zeros(numDirichlet, 1);

for i = 1:numDirichlet
  idx = outerHeart(i, 1);
  cx = geometry(idx, 1);
  cy = geometry(idx, 2);
  
  heartPotential(i, 1) = feInSolution(cx, cy);
end
heartPotential;

% Conductivity on each node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output:   Ndcond = Array of Given Conductivity on each node
% input :   Ndcond = Same of output 
%           heart = Array of Indices of Dirichlet boundary in Coordinate Set

% Conductivity may be given by  \sigma(x,y) = xy or
%                               \sigma(x,y) = 1/(xy) if this is well defined

conductivity = zeros(3, 1, numElements, 4);  
for ei = 1:numElements
  NdinE = connectivities(ei,:);             % [1, 5, 4]
  cxv = geometry(NdinE, 1);                 % [-1 0 -1]'
  cyv = geometry(NdinE, 2);                 % [-1 0  0]'
  
  % \sigma(x,y) =  %%%%%%%% Ex: Function 0;
  if size(Elcond, 1) == 0
    
    condmat =  feInConduct(ei, cxv, cyv);
    for j = 1:3
      conductivity(j, 1, ei, 1) = condmat(j,1,1);
      conductivity(j, 1, ei, 2) = condmat(j,1,2);
      conductivity(j, 1, ei, 3) = condmat(j,1,3);
      conductivity(j, 1, ei, 4) = condmat(j,1,4);
    end
    
  else
    condE = Elcond(ei, 1);
    for j = 1:3
      conductivity(j, 1, ei, 1) = condE;
      conductivity(j, 1, ei, 2) = 0;
      conductivity(j, 1, ei, 3) = 0;
      conductivity(j, 1, ei, 4) = condE;
    end
  end
end

%plot the mesh
if (plotMesh)
  subPlotData(123, numPoints, geometry, connectivities);
end %if

return

function subPlotData(figno, numPoints, geometry, connectivities) 

Xvec = geometry(:,1)';
Yvec = geometry(:,2)';

numElements = size(connectivities, 1);


figure(figno);
triplot(connectivities, Xvec, Yvec, 'Color', [0.9 0.9 0.9]); hold on;

%Plot global nodes
for i = 1:numPoints
  str = num2str(i);
  text(Xvec(1, i), Yvec(1, i), str, 'Color', [0.5 0 0], 'Fontsize', 12);
end

hold on;

return

function [newNodes, newConnectivities] = ...
    subValidateConnectivity(geometry, connectivities)  

numPoints = size(geometry, 1);
numConnect = size(connectivities, 1);

newConnectivities = zeros(numConnect, 3);


ndTgSet = zeros(numPoints, 1);

for ie = 1:numConnect
  for in = 1:3
    node = connectivities(ie, in);
    if ndTgSet(node, 1) == 0
      ndTgSet(node, 1) = 1;
    end
  end
end

newNodes = [];
Map2New = [];

decrf = 0;
newid = 0;
for in = 1:numPoints
  newid = newid + 1;
  
  if ndTgSet(in, 1) == 0
    decrf = decrf + 1;
    Map2New(in, 1) = -1;  % This means this node should be excluded.
  else
    newNodes = [newNodes; geometry(in, :)];
    Map2New(in, 1) = in - decrf;
  end
end


for ie = 1:numConnect
  for j = 1:3
    eid = connectivities(ie, j);
    newConnectivities(ie, j) = Map2New(eid, 1);
  end
end

return
