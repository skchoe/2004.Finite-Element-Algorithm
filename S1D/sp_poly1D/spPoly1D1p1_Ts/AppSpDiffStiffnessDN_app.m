clear

% example commands
left = 0; right = 1;% range of whole domain
% dxvec = [.2, .3, .3, .2]'; %Vector showing length of element
% pvec = [2, 5, 5, 2]'; %Vector showing orders on each element
dxvec = 1*[0.2, 0.2, 0.2, 0.2, 0.2]'; %Vector showing length of element
% pvec = [9, 8, 9, 10, 10]'; %Vector showing orders on each element
pvec = 5*[1,1,1,1,1]';

fFun = 'fKnown';
kFun = 'kKnown';
uFun = 'uSolution';
duFun = 'uDerivat';

vBL = feval(uFun, left);
vBR = feval(duFun, right);

[xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(left, right, dxvec, pvec, vBL, vBR, kFun, fFun);
xsvec = xvec;
yvec_sol = feval(uFun, xsvec);

yvec_sol;
figure(89);
plot(xsvec, yvec_sol, 'r','markersize',5); 
    grid on, title('Ananlytic Solution');

max_error = sp1DCalcErrors(xvec, yvec_apx, yvec_sol)