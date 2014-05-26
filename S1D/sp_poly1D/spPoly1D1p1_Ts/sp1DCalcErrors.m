% max_error = sp1DCalcErrors(yvec_apx, yvec_sol); 

function max_error = sp1DCalcErrors(xvec, yvec_apx, yvec_sol)
s1 = size(yvec_apx, 1);
s2 = size(yvec_sol, 1);

diffvec = yvec_sol-yvec_apx;

figure(93);
plot(xvec, diffvec);
    grid on, title('Error of approximation');

if s1 ~= s2
    Cannot_Compare
else
    max_error = max(abs(diffvec));
end

return