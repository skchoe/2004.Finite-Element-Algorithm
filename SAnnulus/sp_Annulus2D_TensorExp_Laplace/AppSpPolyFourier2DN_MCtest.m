clear


BdyType = 'DN';

dxvec = [.2, .2, .2, .2, .2]'; pvec = [7, 9, 12, 14, 15]';
Nr = size(dxvec,1);
left = .1;
right = left + sum(dxvec);

% example commands
Nth = 32;    % Number of divisions on \theta direction

Ntest = 1;

ranvec = randn(Ntest, 1);
stdrand = std(ranvec)
dataset = [];

G_d = ones(Nth, 1);
G_n = zeros(Nth, 1);

for i=1:Ntest
    
    vdl = 100 + ranvec(i, 1);  vnr = 0;         %Boundary conditions 

    sigleft = ANNEvalUserTensor(left);
    sigright = ANNEvalUserTensor(right);

    [xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec, pvec, Nth, G_d, G_n, sigleft, sigright);

    dataset = [dataset, yvec];
    i
end
% dataset
Nx = size(xvec, 1);

meanvec = zeros(Nx, 1);
stdvec = zeros(Nx, 1);

for j = 1:Ntest
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
% 
% [xvec, yvec_sol] = sp1DEvalExactSol(BdyType, left, right, xvec, vBL, vBR);
% max_error = sp1DCalcErrors(xvec, meanvec, yvec_sol)
