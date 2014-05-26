Readme for spectralab module:

Description:
This program is to verify the spectral element method in 1d
over interval [a, b]. There're programs for spectral method of 1D with 
Dirichlet condition/convergence test and 1D with Dirichlet/Neumann condition
and its convergence test. 

In program, 2 cases of Dirichlet boundary condition also defined
with analytic solutions. 

Function defined:

1. Drawing Approximation and analytic solution with error given order.
   [L, fvector, max_error] = spdiffsineDB(left, right, samples, P, vl, vr, viz)

   	left, right : boundary of domain in 1D
   	vl, vr : boundary condition of u
   	samples : # of ftn values used in drawing a graph
   	P: Order of approximation
   	viz: tag to turn on drawing

   Test examples:

   	spdiffsineDB(-1, 1, 30, 5, 0, 0, 1);
             On [-1,1], with order 5, zero boundary conditions,
             and 30 sample points with result display.

   	spdiffsineDB(0, 2, 30, 5, 1, -1, 1);
             On [0,2], with order 5, boundary conditions u(0)=1,u(1)=-1,
             and 30 sample points with result display.

   	spdiffsineDB(-5, 1, 80, 12, -4, 2, 1);
             On [-50,10], with order 8, u(-50)=-40,u(10)=20,
             and 80 sample points with result display.

2. Drawing graph of errors for each orders-Showing convergence.

   spdiffsineDB_conv(left, right, samples, pvec, vl, vr, viz)

   	left, right : boundary of domain in 1D
   	samples : # of ftn values used in drawing a graph
   	vl, vr : boundary condition of u
   	pvec: Array of Orders
   	viz: tag to turn on drawing

   Test examples:

	pvec = [1:1:12]';

   	spdiffsineDB_conv(-1, 1, 30, pvec, 0, 0, 0);
   	spdiffsineDB_conv(0, 2, 30, pvec, 1, -1, 0);
	spdiffsineDB_conv(-50, 10, 80, pvec, -40, 20, 0);


For routines on Neumann condition, 

Solution in one element
   [L, fvector, max_error] = spdiffsineNB(left, right, samples, P, vl, vr, viz)

   	left, right : boundary of domain in 1D
   	vl : Dirichlet boundary condition of u
	vr : Neumann boundary condition of u
   	samples : # of ftn values used in drawing a graph
   	P: Order of approximation
   	viz: tag to turn on drawing

   spdiffsineNB(-1, 1, 30, 5, 0, 0, 1);
   spdiffsineNB(0, 2, 30, 5, 1, -1, 1);
   spdiffsineNB(-50, 10, 80, 5, -40, 20, 1);
   spdiffsineNB(-0.2, 0.4, 80, 5, 1, 2, 1);

Convergence test

   pvec = [1:2:15]';
   spdiffsineNB_conv(-0.5, 0.3, 30, pvec, -1, 100, 1);
   spdiffsineNB_conv(-1, 1, 30, pvec, 0, 0, 1);
   spdiffsineNB_conv(0, 2, 30, pvec, 1, -1, 1);
   spdiffsineNB_conv(-5, 10.5, 80, pvec, -40, 10, 1);
   spdiffsineNB_conv(-6.5, 6.5, 80, pvec, 0.5, -0.5, 1);
   spdiffsineNB_conv(100.5, 101.4, 40, pvec, 50.5, -1.5, 1);
