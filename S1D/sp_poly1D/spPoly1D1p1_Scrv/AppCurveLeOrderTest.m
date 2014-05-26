%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spAppCurveOrderTest.m Spectral Element Solver on 1D
%
% This is an application module to see how large the order should be 
% to draw a curve having stiff slope in the center of domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

numsampleintervals = 20;
% 
% prob_order = 9;
% [crv_dom, crv_img0,crv_img1,  crv_img2] = spPoissCurveLe(prob_order, numsampleintervals);

NSampleOrders = 3;
Nodr_start = 3;

vdom = zeros(NSampleOrders, 1);
vUleft  = zeros(NSampleOrders, 1);
vUright  = zeros(NSampleOrders, 1);
vDUleft  = zeros(NSampleOrders, 1);
vDUright  = zeros(NSampleOrders, 1);
vD2Uleft  = zeros(NSampleOrders, 1);
vD2Uright  = zeros(NSampleOrders, 1);

for i = 1:NSampleOrders
    prob_order = Nodr_start + 2 * (i-1);
    [crv_dom, crv_img0,crv_img1,  crv_img2] = spPoissCurveLe(prob_order, numsampleintervals);
    
    vdom (i, 1) = prob_order;
    Nright = size(crv_dom, 1);
    vUleft (i, 1) = abs(crv_img0(1, 1));
    vUright (i, 1) = abs(crv_img0(Nright, 1) - 1); 
    vDUleft (i, 1) = abs(crv_img1(1, 1)); 
    vDUright (i, 1) =abs(crv_img1(Nright, 1));
    vD2Uleft (i, 1) = abs(crv_img2(1, 1));
    vD2Uright (i, 1) =abs(crv_img2(Nright, 1));
    
    figure(1970);
    
    subplot(1, 3, 1);
        plot(vdom, vUleft,  'bo-','markersize', 5,'linewidth', 2); hold on;
        plot(vdom, vUright, 'rs--','markersize', 5,'linewidth', 2);
        hold off;
        xlabel('Order of Basis');
        
    subplot(1, 3, 2);    
        plot(vdom, vDUleft,  'bo-','markersize', 5,'linewidth', 2); hold on;
        plot(vdom, vDUright, 'rs--','markersize', 5,'linewidth', 2);
        hold off;
        xlabel('Order of Basis');
        
    subplot(1, 3, 3);
        plot(vdom, vD2Uleft,  'bo-','markersize', 5,'linewidth', 2); hold on;
        plot(vdom, vD2Uright, 'rs--','markersize', 5,'linewidth', 2);
        hold off;
        
%         title('Accuracy of Given RHS function for Boundary Conditions');
        xlabel('Order of Basis');
        
end
    