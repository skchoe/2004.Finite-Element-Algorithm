% function fe2DPoisson(): The FE main solver
%
% Inputs:
%   whichmesh, resolution, bsOdr, degMax
% infoFlag
% Outputs
%   Ag, Fg_1, Fg_2, SolIn, SolBd, BdyDP, BdyDP_2, Ndm,
%   IntNd,BdyN, BdyD, Elm, Hsize



function [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
    feInMeshCons(whichmesh, resolution, BdyDP, Elcond, isFakeSol) 

% Mesh selection and generation/validation of it.
% Save the result into files>>>>>>>>>>>>>>>>>>>

switch whichmesh
  
 case 1
  diV = resolution;
  [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
      feInMeshSqulus(diV);
  Hsize = 1/2^diV;
  
  [gElcond, gBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						gElcond, gBdyDP); 
  
 case 2
  diV = resolution
  leftE = -1;
  rightE = 1;
  [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
      feInMeshSquareDD(leftE, rightE, diV);
  Hsize = (rightE - leftE)/diV
  
  [gElcond, gBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						gElcond, gBdyDP); 
  
 case 3
  diV = resolution
  leftE = 0;
  rightE = 3;
  [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
      feInMeshSquareDN(leftE, rightE, diV);
  Hsize = (rightE - leftE)/diV
  
  [gElcond, gBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						gElcond, gBdyDP); 
  
 case 34    % utahTorso34
  dataFile = 'utahTorso34';
  drawFlag = 0;
  Hsize = 0;
  
  [tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, tElcond] = ...
      feMeshTorsoOrg34(dataFile, Elcond, drawFlag, isFakeSol);
  
  [tElcond, tBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						tElcond, tBdyDP); 
  
  if resolution == 1
    gNdMat = tNdMat; gElm = tElm; gBdyD = tBdyD; gBdyN = tBdyN;
    gIntNd = tIntNd; gBdyDP = tBdyDP; gElcond = tElcond;
    
  elseif resolution == 2
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
        feMeshTorsoBis34(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
			 tElcond, drawFlag, ovFlag); 
    
  elseif resolution == 3
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
        feMeshTorsoTrs34(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
			 tElcond, drawFlag, ovFlag); 
  end %if
  
  
  
 case 37     % utahTorso37
  dataFile = 'utahTorso37_g';
  drawFlag = 0;
  Hsize = 0;
  
  [tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, tElcond] = ...
      feMeshTorsoOrg37(dataFile, Elcond, drawFlag, isFakeSol); 
  
  [tElcond, tBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						tElcond, tBdyDP); 
  
  if resolution == 1
    gNdMat = tNdMat; gElm = tElm; gBdyD = tBdyD; gBdyN = tBdyN; 
    gIntNd; gBdyDP = tBdyDP; gElcond = tElcond;
    
  elseif resolution == 2
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
        feMeshTorsoBis37(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
			 tNdcondF, drawFlag, ovFlag);
    
  elseif resolution == 3
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
        feMeshTorsoTrs37(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
			 tNdcondF, drawFlag, ovFlag);
    
  end %if
  
 case 4
  diV = resolution;
  leftE = -200;
  rightE = 200;
  
  [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
      feInMeshSquareDD(leftE, rightE, diV);
  
  Hsize = (rightE - leftE)/diV;
  
  [gElcond, gBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						gElcond, gBdyDP); 

 case 340    % refinded version of utahTorso34 from CRJ
  dataFile = 'layer70_4';
  drawFlag = 0;
  Hsize = 0;
  
  %should check to see if these have already been processed
  %and avoid processing if so
  [tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, tElcond] = ...
      feMeshTorsoOrg340(dataFile, Elcond, drawFlag, isFakeSol);
  
  [tElcond, tBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						tElcond, tBdyDP); 
   
  if resolution == 1
    gNdMat = tNdMat; gElm = tElm; gBdyD = tBdyD; gBdyN = tBdyN;
    gIntNd = tIntNd; gBdyDP = tBdyDP; gElcond = tElcond;
    
  elseif resolution == 2
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
        feMeshTorsoBis34(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
			 tElcond, drawFlag, ovFlag); 
    
  elseif resolution == 3
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gElcond] = ...
        feMeshTorsoTrs34(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
			 tElcond, drawFlag, ovFlag); 
  end %if

  
end %switch
return

function [tElcond, tBdyDP, ovFlag] = subOverrideInputs(Elcond, BdyDP, ...
						  tElcond, tBdyDP) 

% This step override Conductivity and Dirichlet Bdy Condition from Inputs
ovFlag = 0;
if size(Elcond, 1) ~= 0 % Predefined input
  tElcond = Elcond;     % Override computed gElcond by input Elcond
  ovFlag = 1;
end

if size(BdyDP, 1) ~= 0   % Predefined input
  tBdyDP = BdyDP;         % Override computed Boundary by input
                          % Boundary conditions. 
  ovFlag = 1;
end

return