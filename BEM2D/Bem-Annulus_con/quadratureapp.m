% Testing module for integration using quadrature

% Integrand = u(x) = sin(pi*x), x = [0, 1]

clear
% 
% left = 0;
% right = pi/8;

thi = 4.7124;
thi = 1.1781;
thi = 3*pi/8;
left = -pi/4;  %0.7854;
right = pi/4;  %1.5708;
radius = 1;

Nc = 5;  % Number of tests
St = 10;  % Number of init order

AnsSet = zeros(Nc, 1);
DomSet = zeros(Nc, 1);

for i = 1:Nc
    
    Np = i + St;

    DomSet(i, 1) = Np;
    
    qalpha = 0;
    qbeta = 0;

    % [dgq, dgrq, dglq] = QuadratureDegree(Np)
    [z, w] = JacobiGZW(Np, qalpha, qbeta);
    
    Nzmx = size(z, 1);

    xiz =  (right - left)/2 * z + (left + right)/2 

    J = (right - left) / 2;

    fiz = zeros(Nzmx, 1);
    
    for k = 1:Nzmx
%         fiz(k, 1) = log( sin(xiz(k, 1)) / xiz(k, 1) );  singular alpha1
%         fiz(k, 1) = 2 * log( abs( xiz(k, 1) )); singular alpha2\
%         fiz(k, 1) = log(2 * radius * sin(abs(xiz(k, 1)-thi)/2));
         
%         thx = xiz(k, 1);
%         if thi == thx
%             fiz(k, 1) = log(1);
%         else
%             fiz(k, 1) = log(sin(abs(thx-thi)/2)/(abs(thx-thi)));
%         end

        if xiz(k, 1) == 0
            fiz(k, 1) = 0;
        else
            fiz(k, 1) = log(sin(xiz(k, 1))/xiz(k, 1));
        end
    end
    
%     integral_term = J * dot( fiz, w );
        
%     log_term = (thi-left)*(log(thi-left)-1) + (right-thi)*(log(right-thi)-1);
%     log_term = 0;

    AnsSet(i, 1) =  dot(fiz, w) * J
end

% Ans = 2/pi * ones(Nc, 1);
% plot( DomSet, AnsSet ); hold on
% plot( DomSet, Ans);

