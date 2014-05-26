% FourierDEqn1.m - Approximation of Solution of DE of order 2 Using Fourier Method
%        For complex v, delete "real" commands.

% Since w_hat = i * k * u_hat in differentiation, we use the fact that
% u_hat = -i * 1/k * w_hat, and 1/k = 0 if k==0

  N = 100; h = 2*pi/N; x = h*(1:N)';

% 2-Differentiation of many sample fucntions:
% 1
  f = sin(x); 
  u = - sin(x);
% 2
%   f = (1 - sin(x) - sin(x) .* sin(x)) .* exp(sin(x));
%   u = exp(sin(x));
% 3
%   f = 0; 
%   u = exp(log(x));
% 4
%   f = 1 ./ x; 
%   u = x .* log(x);

  int_ = [0:N/2-1 0 -N/2+1:-1]';  n = size(int_);
  inv_ = zeros(n);  % inv contains reciprocals of 'int' or zeros.
  for i = 1:n
    if int_(i) == 0
        inv_(i) = 0;
    else
        inv_(i) = 1/int_(i);
    end
  end
  
  f_hat = fft(f);
  w_hat = -1 * inv_ .* inv_ .* f_hat;
  w = real(ifft(w_hat));
 
  subplot(3,1,1), plot(x,f,'.-','markersize',13)
    grid on, title('Force Function')

  subplot(3,1,2), plot(x,w,'.-','markersize',13)
    %axis([0 7 -500 500]);
    grid on, title('Fourier Derivative Appx.')

  subplot(3,1,3), plot(x,u,'.-','markersize',13)
    grid on, title('Analytic solution')

  error = norm(w-u,inf)
  text(2.2,1.4,['max error = ' num2str(error)])
