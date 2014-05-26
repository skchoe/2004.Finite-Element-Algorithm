% example commands
left = 0; right = 1;% range of whole domain
dxvec = [1]'; %Vector showing length of element
pvec = [13]'; %Vector showing orders on each element
vdl = 999, vnr = 999;         %Boundary conditions are zero

max_error = sp1DSolvePoisson('DN', left, right, dxvec, pvec, vdl, vnr)
