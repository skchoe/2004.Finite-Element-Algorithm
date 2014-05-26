% example commands
left = -1; right = 1;% range of whole domain
% dxvec = [.2, .3, .3, .2]'; %Vector showing length of element
% pvec = [2, 5, 5, 2]'; %Vector showing orders on each element
dxvec = 2*[0.2, 0.2, 0.2, 0.2, 0.2]'; %Vector showing length of element
pvec = [5, 7, 4, 7, 5]'; %Vector showing orders on each element
vdl = 999; vnr = 999;         %Boundary conditions are zero

BdyType = 'DN';
[xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec, pvec, vdl, vnr);
[xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
max_error = sp1DCalcErrors(yvec_apx, yvec_sol)