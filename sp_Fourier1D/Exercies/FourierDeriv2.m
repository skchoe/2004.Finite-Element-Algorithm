% FourierDeriv2.m - Approximation of Derivative of 2. Using Fourier Method
%        For complex v, delete "real" commands.

% Differentiation of a hat function:
  N = 100; h = 2*pi/N; x = h*(1:N)';

% Differentiation of exp(sin(x)):
%   v = sin(x); 
%   v2prime = - sin(x);

%   v = exp(sin(x)); 
%   v2prime = - sin(x) .* exp(sin(x)) + cos(x).*cos(x).*exp(sin(x));

  v = exp(log(x)); 
  v2prime = 0;
% 
%   v = x.*log(x); 
%   v2prime = 1./x;

  v_hat = fft(v);
  w_hat = -1*[0:N/2-1 0 -N/2+1:-1]'.*[0:N/2-1 0 -N/2+1:-1]' .* v_hat;
  w = real(ifft(w_hat));
 
  subplot(3,1,1), plot(x,v,'.-','markersize',13)
    grid on, title('Original Function')

  subplot(3,1,2), plot(x,w,'.-','markersize',13)
    axis([0 7 -500 500]);
    grid on, title('Spectral 2nd Derivative')

  subplot(3,1,3), plot(x,v2prime,'.-','markersize',13)
    grid on, title('Analytic 2nd Derivative and Error')

  error = norm(w-v2prime,inf)
  text(2.2,1.4,['max error = ' num2str(error)])
