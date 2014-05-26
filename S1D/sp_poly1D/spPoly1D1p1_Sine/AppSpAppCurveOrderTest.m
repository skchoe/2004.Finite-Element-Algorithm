%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spAppCurveOrderTest.m Spectral Element Solver on 1D
%
% This is an application module to see how large the order should be 
% to draw a curve having stiff slope in the center of domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numsampleintervals = 1e+3;

prob_order = 13;
[alphas2, exactsolvec] = spPoissCurve(prob_order, numsampleintervals);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE SOLUTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prob_order == 3
% 
%     Mat = [1, 0, 0, 0;
%             1, 1, 1, 1;
%             0, 1, 0, 0;
%             0, 1, 2, 3];
% 
% prob_order == 5
%     
%     Mat = [1, 0, 0, 0, 0, 0;
%             1, 1, 1, 1, 1, 1;
%             0, 1, 0, 0, 0, 0;
%             0, 1, 2, 3, 4, 5;
%             0, 0, 2, 0, 0, 0;
%             0, 0, 2, 6, 12, 20];
% 
% prob_order == 7    
%     
%     Mat = [1, 0, 0, 0, 0, 0, 0, 0;
%             1, 1, 1, 1, 1, 1, 1, 1;
%             0, 1, 0, 0, 0, 0, 0, 0;
%             0, 1, 2, 3, 4, 5, 6, 7;
%             0, 0, 2, 0, 0, 0, 0, 0;
%             0, 0, 2, 6, 12, 20, 30, 42;
%             0, 0, 0, 6, 0, 0, 0, 0;
%             0, 0, 0, 6, 24, 60, 120, 210];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
