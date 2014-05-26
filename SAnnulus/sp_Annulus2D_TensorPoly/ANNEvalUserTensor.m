


% evaluate tensor function.
function ueval = ANNEvalUserTensor(k, xiz)
    
    dim = size(xiz, 1);
    ueval = ones(dim, 1);
    
    for i = 1:k
        ueval = ueval .* xiz;
    end
    
return