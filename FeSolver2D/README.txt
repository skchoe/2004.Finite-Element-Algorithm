


Readme Jul. 12. 04
- Input argument for controlling quadrature order is updated in all source.
- A execution script "App_FindQdtureOrder.m" is add to find out proper
  quadrature order of numerical integration. It shows quadrature order and
  a graph showing the decrease of sequential error with respect to q-order.

Readme Jul. 09. 04
- Pre/Post processing module is renamed to be distingushed easily.

Readme Jul. 08. 04

This directory is composed of 
- doc: Documentations, such as report of result (Under construction)
- FeSolver2D: Matlab code for 2D Linear solver of Poisson equation.
		App*.m 		- Files for matlab script
		Jacobi*.m - lib.'s for obtaining quadrature points
		fe*.m 		- Finite Element solver for 2D Poisson problem.
			feUser*.m - Files to be defined by a specific problem.


----------------------------------------------------------------------
What has been updated (Jul. 12. 04)?

- Solver using quadratic/cubic basis (global assembly part)

