% utahTorso37_g Info
% Total geometric node = 689
% Nodes which is unreferenced = 0;
% Number of Elements = 1274
% Number of Neumann Boundary : (53) 1-53
% Number of Dirichlet Boundary : (52) 638 - 689
% Number of Interior nodes : (584) 54 - 637




clear all;
close all;

addpath '../FeSolver2D';
drawFlag = 1;
[NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF] = rMeshOrigin('utahTorso37_g', drawFlag);
Hsize = 0.1;

% Mesh refinement from 1 to 4 triangles 
rMeshBisect(NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF);
% % Mesh refinement from 1 to 9 triangles 
rMeshTrisect(NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF);

%----------------------------------------%