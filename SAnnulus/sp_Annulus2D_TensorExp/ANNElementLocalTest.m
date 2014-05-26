%  ANNElementLocalTest
%  Testing script for M1, M2, M3.
order = 5;

a = -1;
b = +1;
%[M1, M2, M3] = ANNElementLocal(a, b, order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since I need z,w's later, I computed JacobiGLZW outside of ANNElementLocal()
%

alpha = 1;
beta = 1;

qalpha = 0; qbeta = 0;  % These are (\alpha, \beta) for quadrature point calculation 
[dgq, dgrq, dglq] = QuadratureDegree(2 * order)
[z, w] = JacobiGZW(dglq, qalpha, qbeta)


Jn = (b-a)/2;
xiz = Jn * z + (b+a)/2

[M1, M2, M3] = ANNElementLocal(order, alpha, beta, z, w, xiz, Jn)
