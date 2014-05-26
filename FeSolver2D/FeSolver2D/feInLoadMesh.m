% function fe2DPoisson(): The FE main solver
%
% Inputs:
%   whichmesh, resolution, bsOdr, degMax
% infoFlag
% Outputs
%   Ag, Fg_1, Fg_2, SolIn, SolBd, BdyDP, BdyDP_2, Ndm,
%   IntNd,BdyN, BdyD, Elm, Hsize



function [Ag, Fg_1, Fg_2, SolIn, SolBd, BdyDP, BdyDP_2, Ndm, ...
	  IntNd, BdyN, BdyD, Elm, Hsize, phiXY, Wpq, JacVec, VecXYm] ...
           = fe2DPoisson(whichmesh, condVec, resolution, bsOdr, ...
			 degMax, isLaplace, infoFlag, itrSize) 

FE_Solver_starts = 1


TimeVec = zeros(8, 1);  % Writing Data/Preprocessing/Quadrature
                        % Points/Basis 
			% Op/LocalMatset/Inversion
			% Generate data and save to files  
COMMENT_WRITING_DATA = 1;
tic;

% Mesh selection and generation/validation of it.
% Save the result into files>>>>>>>>>>>>>>>>>>>
switch whichmesh

  case 1
    diV = resolution;
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
      feInMeshSqulus(diV);
    Hsize = 1/2^diV;

  case 2
    diV = resolution
    leftE = -1;
    rightE = 1;
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
      feInMeshSquareDD(leftE, rightE, diV);
    Hsize = (rightE - leftE)/diV

  case 3
    diV = resolution
    leftE = 0;
    rightE = 3;
    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
      feInMeshSquareDN(leftE, rightE, diV);
    Hsize = (rightE - leftE)/diV

  case 34    % utahTorso34
    dataFile = 'utahTorso34';
    drawFlag = 0;
    Hsize = 0;
    
    [tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, tNdcondF] = ...
      feMeshTorsoOrg34(dataFile, condVec, drawFlag);

    if resolution == 1
      gNdMat = tNdMat; gElm = tElm; gBdyD = tBdyD; gBdyN = tBdyN; gIntNd = tIntNd; gBdyDP = tBdyDP; gNdcondF = tNdcondF;
    
    elseif resolution == 2
      [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
        feMeshTorsoBis34(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
		    tNdcondF, drawFlag); 

    elseif resolution == 3
      [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
        feMeshTorsoTrs34(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, ...
		    tNdcondF, drawFlag); 
    end


  case 37     % utahTorso37
    dataFile = 'utahTorso37_g';
    drawFlag = 0;
    Hsize = 0;
    
    [tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, tNdcondF] = ...
      feMeshTorsoOrg37(dataFile, drawFlag); 
  
    if resolution == 1
      gNdMat = tNdMat; gElm = tElm; gBdyD = tBdyD; gBdyN = tBdyN; gIntNd; gBdyDP = tBdyDP; gNdcondF = tNdondF;
    
    elseif resolution == 2
      [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
        feMeshTorsoBis37(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, tNdcondF);

    elseif resolution == 3
      [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
        feMeshTorsoTrs37(tNdMat, tElm, tBdyD, tBdyN, tIntNd, tBdyDP, tNdcondF);
      
    end
  
  case 4
    diV = resolution;
    leftE = -200;
    rightE = 200;

    [gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF] = ...
      feInMeshSquareDD(leftE, rightE, diV);
  
    Hsize = (rightE - leftE)/diV;
  
end

COMMENT_FINISH_WRITING_DATA = toc;
TimeVec(1,1) = COMMENT_FINISH_WRITING_DATA;

% Load Data, Proc it to Connectivity - Newly arranged by its position
% /Interior/Dboundary/Nboundary
COMMENT_PREPROC_DATA = 2;  
tic;  

[Ndm, Elm, BdyN, BdyD, IntNd, BdyDP, Ndcond, Map2New] = ...
    feDataProc(bsOdr, gNdMat, gElm, gBdyD, gBdyN, gIntNd, gBdyDP, gNdcondF);

%BdyDp = the dirichlet boundary conditions

COMMENT_FINISH_PREPORC_DATA = toc;
TimeVec(2,1) = COMMENT_FINISH_PREPORC_DATA;

COMMENT_PREPARE_STORAGE = 3;
tic;


% For a set of outputs for global assembly
Nnd  = size(Ndm, 1);
NbdD = size(BdyD, 1);
NbdN = size(BdyN, 1);
Nin  = size(IntNd, 1);

Npot = size(BdyDP, 1);
Nco  = size(Ndcond, 1);

Nunkwn  = Nin + NbdN;
SolIn = ones(Nin, 1);
SolBd = ones(NbdN, 1);

% Local Element Stiffness matrices(Upper Diagonal Sequence)
  
Nukwn = NbdN + Nin;                % Nukwn: Number of Unknowns
Ag = sparse(Nukwn, Nukwn);        % Ag: Global Matrix
Fg_1 = sparse(Nukwn, 1);          % Fg_1: First mode (mean) of the
                                  % global RHS  vector
Fg_2 = sparse(Nukwn,1);           % Fg_2: Second mode (variance) of
                                  % the global RHS vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stochastic information for the boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uvec = sparse(Nukwn+NbdD, 1);         % Unknowns+BdyD (mean of the
                                      % boundary conditions)
				      
Uvec_2 = sparse(Nukwn+NbdD,1);        % Second moment of the
                                      % Unknowns +  Boundary Values
Uvec_2(Nukwn+1:Nukwn+NbdD,1) = 0.05;  % Second moment set to
                                      % constant  over all BVs

BdyDP_2 = Uvec_2(Nukwn+1:Nukwn+NbdD,1);


for o = 1:NbdD
  Uvec(Nukwn+o, 1) = BdyDP(o, 1);
end %for o

COMMENT_PREPARESTORAGE = toc;
TimeVec(3,1) = COMMENT_PREPARESTORAGE;


% Local Dimension 'NlocDim'
switch bsOdr
 case 'L'
  NlocDim = 3;
 case 'Q'
  NlocDim = 6;
 case 'C'
  NlocDim = 10;
end

COMMENT_LOCAL_BASIS_COEFF = 4;
tic;

Nel = size(Elm, 1);

% Local Basis Coefficients
coefmat = feLocalBasis(bsOdr);

COMMENT_LOCAL_BASIS_COEFF = toc;
TimeVec(4,1) = COMMENT_LOCAL_BASIS_COEFF;

Nlocmat = size(coefmat, 1);         % 'L':3, 'Q':6, 'C':10
Nloccom = Nlocmat*(Nlocmat+1)/2;    % 'L':6, 'Q':21, 'C':55 
				    % Number of upper diagonal local elements
				    
% LocStMat : rowxcolumn - Array of Nlocmat x Nlocmat local stiffness matrix
%            3rd dim - Number of local elements 
% LocSoMat : rowx1 - Array of Nlocmat x 1 local source matrix 
LocStMat = zeros(Nlocmat, Nlocmat);
LocSoMat = zeros(Nlocmat, 1);


COMMENT_QUADRATUREPOINTS = 5;
tic;

% Quadrature values
[zwTableL, zwTableR] = subFindQuadratures(degMax);
NzL = size(zwTableL, 1);
NzR = size(zwTableR, 1);

COMMENT_END_QPOINTS = toc;
TimeVec(5,1) = COMMENT_END_QPOINTS;

COMMENT_COMMON_OP_OVER_ELTS = 6;    
tic;

% Compute quadrature in canonical space(xi1, xi2)
VecXY = zeros(NzL, NzR, 2);
VecXY = subXformToRefSpace(NzL, NzR, zwTableL, zwTableR, VecXY);


% Basis table evaluated at each element of zwTableL and zwTableR
phiXY = zeros(NzL, NzR, NlocDim);
phidX = zeros(NzL, NzR, NlocDim);
phidY = zeros(NzL, NzR, NlocDim);

phXYW = zeros(NzL, NzR, NlocDim);
Wpq = zwTableL(:,2) * zwTableR(:,2)';

% Evaluate Basis functions on quadrature points
[phiXY, phidX, phidY, phXYW] = subBasisPreEvaluation ...
    (bsOdr, NlocDim, NzL, NzR, coefmat, Wpq, VecXY, phiXY, phidX, ...
     phidY, phXYW); 

% Preprocessing for local stiffness matrix components(products
% between basis function values) 
PhiS11 = zeros(NzL, NzR, NlocDim, NlocDim);
PhiS12 = zeros(NzL, NzR, NlocDim, NlocDim);
PhiS21 = zeros(NzL, NzR, NlocDim, NlocDim); % Not computed for
                                            % symmetry of matrix 
PhiS22 = zeros(NzL, NzR, NlocDim, NlocDim);

[PhiS11, PhiS12, PhiS21, PhiS22] = subBasisPreProduct ...
    (NlocDim, NzL, NzR, phidX, phidY, PhiS11, PhiS12, PhiS21, PhiS22);

COMMENT_END_COMMON_OP_OVER_ELTS = toc;
TimeVec(6,1) = COMMENT_END_COMMON_OP_OVER_ELTS;


COMMENT_LOCAL_MAT = 7;
tic;

%%%%% PER ELEMENT OPERATIONS %%%%%

% Quadrature points in physical space: a matrix for all elements
VecXYm = zeros(NzL, NzR, Nel, 2);
JacVec = zeros(Nel, 1);

for ie = 1:Nel % For each element,
  
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie==1 
    tic;    % differential, jacobian
  end
  
  
  [jcb, dxxi] = subLocalXform(Ndm, Elm(ie, :)); % Note that jcb(jacobian) 
						% is not an absolute value
						% dxxi = d(xi_1,
                                                % xi_2)/d(x_1,
                                                % x_2)^T  
						% = 1/jcb *
                                                % [dx2/dxi2
                                                % -dx2/dxi1;-dx1/dxi2
                                                % dx1/dxi1]; 
						% inside of subLocalXform()
						JacVec(ie, 1) = abs(jcb);
						
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie == 1
    T1_jacob_difl = toc; % differential, jacobian
    tic;                 % precomputing of products of differentials
  end
  
  
  % Element specific preparation for evaluating integrand
  % 1. Differential for (xi1,xi2) w.r.t. (x1, x2)
  RefD = dxxi * transpose(dxxi);
  
  DifMat2 = 1/jcb^2* [RefD(4,4) -RefD(4,3) -RefD(3,4)  RefD(3,3);
		      -RefD(2,4)  RefD(2,3)  RefD(1,4) -RefD(1,3);
		      -RefD(4,2)  RefD(4,1)  RefD(3,2) -RefD(3,1);
		      RefD(2,2) -RefD(2,1) -RefD(1,2)  RefD(1,1)];
  
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie == 1
    T2_multdifl = toc;   % precomputing of products of differentials
    tic;                % coordinate transform for RHS integrand evaluation
  end
  
  % 2. Quadrature points from reference coord. to physical coord.
  coord = zeros(3, 2);    % (x,y)Coordinate of 3 nodes of the element
%   condc = zeros(3, 4);    % 2x2 = 4x1 Conductivity for 
                           % 3-coordinates on this element
  
  for s = 1:3
    coord(s, 1) = Ndm(Elm(ie, s), 2);   
    coord(s, 2) = Ndm(Elm(ie, s), 3);
%     condc(s, :) = Ndcond(Elm(ie, s), :);
  end
  
  % Quadrature points in physical space:
  % (VecX, VecY) -> (VecXe,VecYe) \in physic. sp.
  % Interpolation factors (fac1, fac2, fac3) from (xi1, xi2) to (x1, x2)
  VecXYe = zeros(NzL, NzR, 2);
  % Conductivity on each quadrature points: 4 components of 2x2 matrices
  Cond4 = zeros(NzL, NzR, 4);
  
  for iq = 1:NzR
    tmX = VecXY(:,iq, 1);
    tmY = VecXY(:,iq, 2);
    fac1 = (-tmX - tmY)/2;
    fac2 = ( 1 + tmX)/2;
    fac3 = ( 1 + tmY)/2;
    
    % Coordinate Xform to (x_1, x_2) <-Need in mass matrix, not in
    % stiffness matrix 
    for ci = 1:2    % x, y coordinates
      VecXYe(:, iq, ci) = coord(1,ci) * fac1 + coord(2,ci) * fac2 + ...
	  coord(3,ci) * fac3; 
    end
    
    % Conductivity should obtained from 1) exact function type or
    % result of interpolation from known values on nodes 
    %           for cj = 1:4    % 4 components showing tensor on
    %                           % (NzL, NzR) nodes 
    %               Cond4(:, iq, cj) = condc(1,cj) * fac1 +
    %               condc(2,cj) * fac2 + condc(3,cj) * fac3; 
    %           end            
    Cond4 = feInConduct(ie, VecXYe(:,:,1), VecXYe(:,:,2));
    
  end
  
  % Collect quadrature points in physical space into VecXYm (to get
  % them for all elements) 
  VecXYm(:, :, ie, 1) = VecXYe(:, :, 1);
  VecXYm(:, :, ie, 2) = VecXYe(:, :, 2);
  
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie == 1
    T3_qptevalvec = toc;     % coordinate transform for RHS
                             % integrand evaluation 
    tic;                     % evaluation of source on each
                             % quadrature coordinate 
  end%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  % Evaluate Source over (x1, x2)
  % In case of solver's being Laplace solver, we set source 
  % term to be zero all over the quadrature points;
  if isLaplace == 1
    eF = zeros(NzL, NzR);
  else
    %eF = feval(@feInSource, ie, VecXYe(:,:,1), VecXYe(:,:,2));
    eF = feInSource(ie, VecXYe(:,:,1), VecXYe(:,:,2));
  end
  
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie == 1
    T4_sourceeval = toc;     % evaluation of source on each
                             % quadrature coordinate 
    tic;                     % local stiffness matrix, RSH vector evaluation
  end%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  for j = 1:Nlocmat
       
    if isLaplace == 1
      LocSoMat(j, 1) = 0;
    else
      LocSoMat(j, 1) = sum(sum( eF .* phXYW(:,:,j) * abs(jcb)));
    end
    
    for k = j:Nlocmat % For each (j, k) component of local matrix
      
      locIntMat = zeros(4,4);
      
   for p = 1:NzL
	for q = 1:NzR
	  PhiSvec = [PhiS11(p,q,j,k) PhiS12(p,q,j,k) PhiS21(p,q,j,k) ...
		     PhiS22(p,q,j,k)]; % PhiS12(:,:,j,k) == PhiS21(:,:,k,j) 
	  PhiSvec = Wpq(p,q) * PhiSvec;
	  Cdvec = [Cond4(p,q,1);Cond4(p,q,2);Cond4(p,q,3);Cond4(p,q,4)];
	  locPQ = Cdvec * PhiSvec;
	  locIntMat = locIntMat + locPQ;
	end
      end
      
      LocStMat(j, k) = sum(sum(DifMat2 .* locIntMat)) * abs(jcb);
      
      if k ~= j
	LocStMat(k, j) = LocStMat(j, k);
      end
      
    end
  end
  
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie == 1
    T5_localmat = toc;  % local stiffness matrix, RSH vector evaluation
    tic;                 % unit-global assemblly(per element)
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % GLOBAL ASSEMBLY
  
  %       ie
  %       globalnode = Elm(ie, :)
  
  for row = 1:NlocDim
    
    grow = Elm(ie, row);
    
    
    if grow <= Nukwn
      
      Fg_1(grow, 1) = Fg_1(grow, 1) + LocSoMat(row, 1);
      Fg_2(grow, 1) = 0;
      
      for col = 1:NlocDim
	
	gcol = Elm(ie, col);
	LocElt = LocStMat(row, col);
	
	if gcol <= Nukwn
	  Ag(grow, gcol) = Ag(grow, gcol) + LocElt;
	else
	  Fg_1(grow, 1) = Fg_1(grow, 1) - LocElt * (Uvec(gcol, 1));
	  Fg_2(grow, 1) = Fg_2(grow, 1) - LocElt * Uvec_2(gcol,1);
	end
	
      end
      
    end
    
  end
  
  
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie == 1
    T6_globalassm = toc;     % unit-global assemblly(per element)
  end
  
  
  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if ie == 1
    
    if infoFlag == 1
      elTimeVecInfo = ['1_jacobian';'2_difprdct';'3_qptevalv'; ...
		       '4_sourcevl';'5_locmtset';'6_insrglmt'] 
      elTimeVec = [T1_jacob_difl;T2_multdifl;T3_qptevalvec; ...
		   T4_sourceeval;T5_localmat;T6_globalassm] 
      %           figure(55);
      %           stem(eTimeVec);
      TotalLocCompT = sum(elTimeVec)
    end
  end
  
end

COMMENT_LOCAL_MAT = toc;
TimeVec(7,1) = COMMENT_LOCAL_MAT;

COMMENT_INVERSION = 8;
tic;

% detA = det(Ag);
% if detA ~=0
  Nusol = sparse(Nunkwn, 1);
  Nusol = pcg(Ag, Fg_1, 10e-12, itrSize);
%   Nusol = Ag \ Fg_1;  

  for i = 1:Nin
    SolIn(i, 1) = Nusol(i, 1);
  end
  
  for j = 1:NbdN
    SolBd(j, 1) = Nusol(j+Nin, 1);
  end
% else
%   DETERMINANT_ZERO = 0
% end

COMMENT_INVERSION = toc;
TimeVec(8,1) = COMMENT_INVERSION;

if infoFlag == 1 
  TimeVecInfo = ['WRTINGDATA';'PREPROCDAT';'PREPARSTOR';'LOCBASISCF';'QDRTPOINTS';'ED_COMMELT';'LOCAL_MATX';'INVER_SION']
  TimeVec
end

return


function [zwTableL, zwTableR] = subFindQuadratures(degMax)
    % Finding a table of quadrature points/weights for each element%%%%%%%%%
    % And Basis value table for each element %%%%%%%%%%%%%%%%%%%%
    qalpha1 = 0; qbeta1 = 0;  % These are (\alpha, \beta) for
                              % quadrature point calculation  
    qalpha2 = 1; qbeta2 = 0;  % These are (\alpha, \beta) for
                              % quadrature point calculation  

    [dgq, dgrq, dglq] = QuadratureDegree(degMax);   % find number
                                                    % of q-pts with
                                                    % respect to
                                                    % degree
                                                    % 'degMax' of
                                                    % integrand. 
    [zL, wL] = JacobiGLZW(dglq, qalpha1, qbeta1);
    [zR, wR] = JacobiGRZW(dgrq, qalpha2, qbeta2);
    wR = wR/2;

    zwTableL = zeros(size(zL, 1), 2);
    zwTalbeR = zeros(size(zR, 1), 2);
    zwTableL(:, 1) = zL;
    zwTableL(:, 2) = wL;
    zwTableR(:, 1) = zR;
    zwTableR(:, 2) = wR;

return



% Jacobian d(x_1,x_2)/d(\xi_1, \xi_2) between (x_12), (\xi_12) coordinates
% Based on equation (4.36)
function [jcb, dxxi] = subLocalXform(Ndm, Ndids)

    % Coordinates that forms ith local element
    id1 = Ndids(1, 1);
    id2 = Ndids(1, 2);
    id3 = Ndids(1, 3);
    
    % x1, y1, x2, y2, x3, y3
    x1 = Ndm(id1, 2);   y1 = Ndm(id1, 3);
    x2 = Ndm(id2, 2);   y2 = Ndm(id2, 3);
    x3 = Ndm(id3, 2);   y3 = Ndm(id3, 3);

    dxxi = zeros(4, 1);
    dxxi(1,1) = -.5 * x1 + .5 * x2; 
    dxxi(2,1) = -.5 * x1 + .5 * x3;
    dxxi(3,1) = -.5 * y1 + .5 * y2;
    dxxi(4,1) = -.5 * y1 + .5 * y3;
    
    jcb = dxxi(1,1) * dxxi(4,1) - dxxi(2,1) * dxxi(3,1);
    
return

% coord = (2, 1): (x_1; x_2) for each coordinate
% bcth = coefficient set of local basis function
% dvec  = (2 ,1): gradient
function vec = subEvalBasis(bsOdr, bcth, coord)

    cx = coord(1, 1);
    cy = coord(2, 1);

    switch bsOdr
      case 'L'
        e = [1; cx; cy];
      case 'Q'
        e = [1; cx; cy; cx*cx; cx*cy; cy*cy];
      case 'C'
        e = [1; cx; cy; cx*cx; cx*cy; cy*cy; cx*cx*cx; cx*cx*cy; ...
	     cx*cy*cy; cy*cy*cy]; 
    end

    vec = dot(bcth, e);

return


% coord = (2, 1): (x_1; x_2) for each coordinate
% bcth = coefficient set of local basis function
% dxxi = partial derivative of x for \xi->considered next routine, not here.
% dvec  = (2 ,1): gradient
function dvec = subEvalBasisDeriv(bsOdr, bcth, coord)

    cx = coord(1, 1);
    cy = coord(2, 1);

    switch bsOdr
      case 'L'
        e1 = [0; 1; 0];
        e2 = [0; 0; 1];
      case 'Q'
        e1 = [0; 1; 0; 2 * cx; cy; 0];
        e2 = [0; 0; 1; 0; cx; 2 * cy];
      case 'C'
        e1 = [0; 1; 0; 2 * cx; cy; 0; 3 * cx^2; 2 * cx * cy; cy^2; 0];
        e2 = [0; 0; 1; 0; cx; 2 * cy; 0; cx^2; 2 * cx * cy; 3 * cy^2];
    end

    dvec(1, 1) = dot(bcth, e1);
    dvec(2, 1) = dot(bcth, e2);
   
return

% Transform quadrature points to reference triangluar space
function  VecXY = subXformToRefSpace(NzL, NzR, zwTableL, zwTableR, VecXY)

  for p = 1:NzL
      xiL = zwTableL(p, 1);
      for q = 1:NzR
          xiR = zwTableR(q, 1);
          % Coordinate Xform to (\xi1, \xi2)
          VecXY(p, q, 1) = (1+xiL).*(1-xiR)./2 - 1;
          VecXY(p, q, 2) = xiR;
      end
  end

return


function [phiXY, phidX, phidY, phXYW] ...
= subBasisPreEvaluation(bsOdr, NlocDim, NzL, NzR, coefmat, Wpq, ...
			VecXY, phiXY, phidX, phidY, phXYW) 
 
  for k = 1:NlocDim

    bcJ = coefmat(:, k);
    
    for q = 1:NzR
      for p = 1:NzL
	xiL = VecXY(p, q, 1);
	xiR = VecXY(p, q, 2);
	% Evaluate local basis on each node in canonical coordinate
	phiXY(p, q, k) = subEvalBasis(bsOdr, bcJ, [xiL;xiR]);
	% Evaluate derivatives of local basis on each node in
	% canonical coordinate 
	dxi12 = subEvalBasisDeriv(bsOdr, bcJ, [xiL;xiR]);
	phidX(p, q, k) = dxi12(1, 1);
	phidY(p, q, k) = dxi12(2, 1);
      end
    end
    
    % RHS precomputation
    phXYW(:,:,k) = phiXY(:,:,k).*Wpq;
    
  end
  
return


function [PhiS11, PhiS12, PhiS21, PhiS22] = subBasisPreProduct(NlocDim, ...
						  NzL, NzR, phidX, ...
						  phidY, PhiS11, ...
						  PhiS12, PhiS21, PhiS22)   

  for i = 1:NlocDim
    for j = 1:NlocDim
      for p = 1:NzL
	for q = 1:NzR
	  PhiS11(p, q, i, j) = phidX(p, q, i) * phidX(p, q, j);
	  PhiS22(p, q, i, j) = phidY(p, q, i) * phidY(p, q, j);
	  PhiS12(p, q, i, j) = phidX(p, q, i) * phidY(p, q, j);
	  PhiS21(p, q, i, j) = phidY(p, q, i) * phidX(p, q, j); %
                                                                % Not
                                                                % computed
                                                                % for
                                                                % symmetry
                                                                % of
                                                                % matrix  
	end
      end                  
    end
  end
  
return

