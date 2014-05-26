##########################################################################
README of Monte Carlo testing of 1D Spectral Poly Method
Feb. 11. 04
##########################################################################
This is test for eqn.

d/dx C d/dx u(x) = sin(pi*x), x \in [-1,1] 

This has exact solution U.

Boundary conditions
u(-1) = 100, d/dx u(1) = 0

Num. MC iteration 
1000

Constant tensor
C = 5

Sample info
1) N(0, 1)
2) Standard deviation: 0.9526

Output info

1) MeanVarCrv.eps : 
	-Mean of solutions(Blue line)
	-Upper bdy of standard deviation from Mean(Red line)
	-Lower bdy of standard deviation from Mean(Magenta line)

2) ErrToExactSol.eps:
        -A graph of curve of difference 
 	 between Mean of solutions and Exact solution.

3) Total Error
	-Max norm of 2)
 	-1.3157

4) Comparison in Standard Deviation between Boundary samples and Solutions in
each 
