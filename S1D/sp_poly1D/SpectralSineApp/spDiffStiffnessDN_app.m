% example commands
left = 0; right = 1;% range of whole domain
dxvec = [1]'; %Vector showing length of element
pvec = [7]'; %Vector showing orders on each element
vdl = 999, vnr = 999;         %Boundary conditions are zero

max_error = spDiffStiffnessDN(left, right, dxvec, pvec, vdl, vnr, numsampleintervals)
