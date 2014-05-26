clear

% example commands
left = -1; right = 1;% range of whole domain
% dxvec = [.2, .3, .3, .2]'; %Vector showing length of element
% pvec = [2, 5, 5, 2]'; %Vector showing orders on each element
dxvec = 2*[0.2, 0.2, 0.2, 0.2, 0.2]'; %Vector showing length of element
pvec = [9, 8, 9, 10, 10]'; %Vector showing orders on each element
% pvec = [1, 2, 1, 2, 1]'; %Vector showing orders on each element

BdyType = 'DN';


Ntest = 1000;

ranvec = randn(Ntest, 1);
stdrand = std(ranvec)
dataset = [];

for i=1:Ntest
    
    vdl = 100 + ranvec(i, 1);  vnr = 0;         %Boundary conditions 

    sigleft = spEvalUserTensor(left);
    sigright = spEvalUserTensor(right);

    [xvec, yvec_apx, vBL, vBR] = sp1DSolvePoisson(BdyType, left, right, dxvec, pvec, vdl, vnr, sigleft, sigright);
    % [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
    % max_error = sp1DCalcErrors(xvec, yvec_apx, yvec_sol);
    
    dataset = [dataset, yvec_apx];
    i
end
% dataset
Nx = size(xvec, 1);

meanvec = zeros(Nx, 1);
stdvec = zeros(Nx, 1);

for j = 1:Nx
    X = dataset(j, :);
    meanvec(j, 1) =  mean(X);
    stdvec(j, 1) = std(X);
end


figure(43);
plot(xvec, meanvec, 'bo-','markersize', 5,'linewidth', 2); hold on;
plot(xvec, meanvec + stdvec/2, 'rs--','markersize', 3,'linewidth', 1); hold on;
plot(xvec, meanvec - stdvec/2, 'ms--','markersize', 3,'linewidth', 1);
    grid on;
    xlabel('x-Domain (Gauss Jacobi Points)');
    ylabel('Distribution of the Solution');
legend('Mean Values','Mean + Standard Deviation/2', 'Mean - Standard Deviation/2');

[xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
max_error = sp1DCalcErrors(xvec, meanvec, yvec_sol)
