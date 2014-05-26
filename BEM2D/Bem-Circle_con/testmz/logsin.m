% input is row vector when this function is used as an argument of quad();
% output is also row vector form - quad()

% input is column vector when this function is used as an argument of feval();
% output is column vector - feval();

function y=logsin(x)

    X = x';
    N = size(X, 1);
    Y = [];
    for i=1:N
        elt = X(i, 1);
        
        if elt == 0
            frac = 1;
        else
            frac = sin(abs(elt))./abs(elt);
        end
        
        tmp = log(frac);
       
        Y = [Y;tmp];
    end
    
    y = Y';
    
return