clear


BdyType = 'DN';

dxvec = [.2, .2, .2, .2, .2]'; pvec = [7, 9, 12, 14, 15]';

Nr = size(dxvec,1);
left = .1;
right = left + sum(dxvec);

% example commands
Nth = 32;    % Number of divisions on \theta direction

Ntest = 1000;

ranvec = randn(Ntest, 1);
stdrand = std(ranvec)

deriv_l = 0;
deriv_r = 1;

G_d = ones(Nth, 1);
G_n = zeros(Nth, 1);

ueval = zeros(Nth, 2);
ueval = ANNEvalUserSolution(left, right, Nth, deriv_l, deriv_r, ueval);

for i=1:Ntest
    
    vdl = 100 + ranvec(i, 1);  vnr = 0;         %Boundary conditions 

    sigleft = ANNEvalUserTensor(left);
    sigright = ANNEvalUserTensor(right);

    [xvec, yvec, zMat_apx] = ANN2DSolvePoisson(BdyType, left, right, dxvec, pvec, Nth, ueval(:, 1), ueval(:, 2), sigleft, sigright);

    if i==1
        dataset = zeros(size(xvec, 1), size(yvec,  1), Ntest)
    end
    SN = size(zMat_apx)
    dataset(:, :, i) = zMat_apx';
    dataset
    i
end
% dataset
Nx = size(xvec, 1);
Ny = size(yvec, 1);

mean_mat = zeros(Nx, Ny);
std_mat = zeros(Nx, Ny);

for j = 1:Nx
    for k = 1:Ny
        X = dataset(j, k, :);
        mean_mat(j, k) =  mean(X);
        std_mat(j, k) = std(X);
    end
end

figure(43);

mMat = real(mean_mat);
% s1 = size(xvec)
% s2 = size(yvec)
% s3 = size(mMat)
% mesh(yvec, xvec, mMat);

mesh(yvec, xvec, mMat,'markersize', 5,'linewidth', 2); hold on;
% mesh(yvec, xvec, mMat, 'bo-','markersize', 5,'linewidth', 2); hold on;
mesh(yvec, xvec, mMat + 10 * std_mat/2, 'markersize', 3,'linewidth', 1); hold on;
mesh(yvec, xvec, mMat - 10 * std_mat/2, 'markersize', 3,'linewidth', 1);
    grid on;
    xlabel('y-Domain (Gauss Jacobi Points)');
    ylabel('x-Domain (Gauss Jacobi Points)');
 legend('Mean Values','Mean + Standard Deviation/2', 'Mean - Standard Deviation/2');
