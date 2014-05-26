% Script M-file to test function spEvalMultiSpan1D


xvec = [-1;0;1];
eleftvec = -1;
erightvec = 1;
alpha = 1; 
beta = 1;
ug_hat = [1, 0, 1];
pvec = [1,2,3];
map = [1,2,3;1,2,3];
ethvec = [1,2,3];
%
yvec = spEvalMultiSpan1D(xvec, eleftvec, erightvec, alpha, beta, ug_hat, pvec, map, ethvec);
yvec
%
