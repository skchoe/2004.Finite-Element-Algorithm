clear

% example commands
% left = 0; right = 1;% range of whole domain
% dxvec = [.25, .25, .25, .25]'; %Vector showing length of element
% pvec = [2, 2, 2, 2]'; %Vector showing orders on each element
% vdl = 999; vdr = 999;         %Boundary conditions are zero
% numsampleintervals = 1000;
% left = 0; right = 0.00001; dxvec = [0.0000025, 0.0000025, 0.0000025, 0.0000025]'; pvec = [10, 10, 10, 10]'; vdl = 999; vdr = 999; numsampleintervals = 1000;
% left = 0; right = 1e-5; dxvec = [2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6]'; pvec = [10, 10, 10, 10]'; vdl = 999; vdr = 999; numsampleintervals = 1000;
%left = 0; right = 1; dxvec = [2.5e-1, 2.5e-1, 2.5e-1, 2.5e-1]'; pvec = [21, 10, 10, 10]'; vdl = 999; vdr = 999; numsampleintervals = 1000;
p=9;
left = 0; right = 1; dxvec = [.3, .4, .3]'; pvec = p*[1, 1, 1]'; vdl = 999; vdr = 999;

max_error = sp1DSolvePoisson('DD', left, right, dxvec, pvec, vdl, vdr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FTN1  left = 0.0; right = 1.0; dxvec = [1]'; pvec = [15]'; vdl = 999; vdr = 999; numsampleintervals = 1000; ->Error 1e-16
%FTN3  left = 0.0; right = 1.0; dxvec = [.5, .5]'; pvec = [2,2]'; vdl = 999; vdr = 999; numsampleintervals = 1000; ->Error 0.0475
%FTN3  left = 0; right = 1; dxvec = [.25, .25, .25, .25]'; pvec = [2, 2, 2, 2]'; vdl = 999; vdr = 999; numsampleintervals = 1000; ->Error 0.0098
%FTN3  left = 0; right = 1e-5; dxvec = [2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6]'; pvec = [10, 10, 10, 10]'; vdl = 999; vdr = 999; numsampleintervals = 1000;->Error  1.3452e-043
%left = 0; right = 1; dxvec = [1]'; pvec = [15]'; vdl = 999; vdr = 999; numsampleintervals = 10000; error(p=15):5.173139694392148e-011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
