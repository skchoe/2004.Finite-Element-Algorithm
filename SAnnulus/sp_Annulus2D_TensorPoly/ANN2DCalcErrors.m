% max_error = ANN2DCalcErrors(zMat_apx, zMat_sol)

function max_error = ANN2DCalcErrors(zMat_apx, zMat_sol)
    s1 = size(zMat_apx);
    s2 = size(zMat_sol);

    if s1 ~= s2
    else
        max_error_mat = max(abs(zMat_sol - zMat_apx));
    end
    
    max_error = max(max_error_mat)
return