% input is row vector when this function is used as an argument of quad();
% output is also row vector form - quad()

% input is column vector when this function is used as an argument of feval();
% output is column vector - feval();

function y = fn_logsin(x)

    X = x'; % column vector conversion
    center = pi/16;
    N = size(X, 1);
    Y = [];

    for i=1:N
        elt = abs(X(i, 1) - center)/2;

        if elt == 0
            frac = 1;
        else
            frac = sin(elt)./elt;
        end
        
        tmp = log(frac);
       
        Y = [Y;tmp];
    end
    
    y = Y';
    
return