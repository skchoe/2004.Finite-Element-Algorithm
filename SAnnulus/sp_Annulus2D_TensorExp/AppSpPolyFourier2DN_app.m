clear

p = 3;
% dxvec = [1]'; pvec = [p]';    %p*[1]';
% 
% dxvec = [.5, .5]'; pvec = [p, p+1]';    %p*[1]';
% dxvec = [.5, .5]'; pvec = p*[1, 1]';
% dxvec = [.1, .1, .1, .1, .1, .1, .1, .1, .1, .1]'; pvec = p*[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]';
% dxvec = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]'; pvec = [7, 7, 14, 15, 13, 13, 15, 14, 7, 7]';
dxvec = [.2, .2, .2, .2, .2]'; pvec = [7, 9, 12, 14, 15]';
Nr = size(dxvec,1);
left = .1;
right = left + sum(dxvec);

% example commands
Nth = 32;    % Number of divisions on \theta direction

% Polynomial/Fourier Method Call
% (xvec, yvec): Grid of domain values;
% zMat_apx, zMat_sol: Approximated / Analytic solutions;
% vBL, vBR: final form of vectors showing boundary values;

BdyType = 'DN';

% Boundary vectors preprocessing
deriv_l = 0;    % degree of derivative for left  bdy condition
deriv_r = 1;    % degree of derivative for right bdy condition

ueval = zeros(Nth, 2);
ueval = ANNEvalUserSolution(left, right, Nth, deriv_l, deriv_r, ueval);

sigleft = ANNEvalUserTensor(left);
sigright = ANNEvalUserTensor(right);

[xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec, pvec, Nth, ueval(:, 1), ueval(:, 2), sigleft, sigright);

[xvec, yvec, zMat_sol] = ANN2DEvalExactSol(left, right, xvec, yvec);

max_error = ANN2DCalcErrors(zMat_apx, zMat_sol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FTN1  left = 0.0; right = 1.0; dxvec = [1]'; pvec = [15]'; vdl = 999; vdr = 999; numsampleintervals = 1000; ->Error 1e-16
%FTN3  left = 0.0; right = 1.0; dxvec = [.5, .5]'; pvec = [2,2]'; vdl = 999; vdr = 999; numsampleintervals = 1000; ->Error 0.0475
%FTN3  left = 0; right = 1; dxvec = [.25, .25, .25, .25]'; pvec = [2, 2, 2, 2]'; vdl = 999; vdr = 999; numsampleintervals = 1000; ->Error 0.0098
%FTN3  left = 0; right = 1e-5; dxvec = [2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6]'; pvec = [10, 10, 10, 10]'; vdl = 999; vdr = 999; numsampleintervals = 1000;->Error  1.3452e-043
%left = 0; right = 1; dxvec = [1]'; pvec = [15]'; vdl = 999; vdr = 999; numsampleintervals = 10000; error(p=15):5.173139694392148e-011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%