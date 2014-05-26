


% evaluate tensor function.

%   Ex 3: Sine Curve
%function ueval = ANNEvalUserTensor(k, xiz)
    
%    dim = size(xiz, 1);
%    ueval = ones(dim, 1);
    
%    for i = 1:k
%        ueval = ueval .* xiz;
%    end
%    
%return

%   Ex 4: S-Curve
function ueval = ANNEvalUserTensor(xiz)

onevec = ones(size(xiz, 1), 1);
ueval = exp( - (xiz - onevec) .* (xiz - onevec) );


return
