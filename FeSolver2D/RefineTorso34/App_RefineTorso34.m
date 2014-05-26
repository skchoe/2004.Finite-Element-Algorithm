% utahTorso34 Info
% Total geometric node = 659 (But 659th node is never referenced)
% Nodes which is unreferenced = 1;
% Number of Elements = 1203
% Number of Neumann Boundary : (53) 1-51, 580, 650
% Number of Dirichlet Boundary : (60) 198 - 257
% Number of Interior nodes : (545) {52, ...,197,258,...579,581,...,649,651,...658} 

clear all;
close all;

% addpath '../FeSolver2D';
drawFlag = 0;
[NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF] = rMeshOrigin34('utahTorso34', drawFlag);
Hsize = 0.1;

org_Nnd = size(NdMat, 1)
org_Nelt = size(Elm, 1)
org_NbdD = size(BdyD, 1)
org_NbdN = size(BdyN, 1)
org_NintN = size(IntNd, 1)

% Mesh refinement from 1 to 4 triangles 
drawFlag = 0;
bis_spec = rMeshBisect34(NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF, drawFlag)


% % % Mesh refinement from 1 to 9 triangles 
drawFlag = 0;
tri_spec = rMeshTrisect34(NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF, drawFlag)

%----------------------------------------%