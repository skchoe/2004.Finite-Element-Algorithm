Np = 22;
qalpha = 0;
qbeta = 0;
[z, w] = JacobiGLZW(Np, qalpha, qbeta);
Nzmx = size(z, 1);


% int{-1,1} ln|x| dx
z
integral1 = dot(log(abs(z)), w)

integral2 = -2