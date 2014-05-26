


% evaluate tensor function.
function ueval = ANNEvalUserTensor(k, xiz)
    
    %   Note that xiz != 0
    ueval = k * xiz^(k-1);

return