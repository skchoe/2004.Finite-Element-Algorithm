% ch3_1.m - Exercies of Chapter 3

No_ = 1
% Set up grid and differentiation matrix:
  N = 2
  h = 2*pi/N; 
  x = h*(1:N)';
  column1 = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
  D1 = toeplitz(column1,column1([1 N:-1:2]))
  E1 = D1*D1
  F1 = E1*D1
  column2 = [-pi^2/(3*h^2)-1/6 .5*(-1).^(1:N-1)./(sin((1:N-1)*h/2))^2]';
  D2 = toeplitz(column2,column2([1 N:-1:2]))
  Are_D2s_equal = isequal(E1, D2)
  
  N = 4
  h = 2*pi/N; 
  x = h*(1:N)';
  column1 = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
  D1 = toeplitz(column1,column1([1 N:-1:2]))
  E1 = D1*D1
  F1 = E1*D1
  column2 = [-pi^2/(3*h^2)-1/6 .5*(-1).^(1:N-1)./(sin((1:N-1)*h/2)).^2]'
  D2 = toeplitz(column2,column2([1 N:-1:2]))
  Are_D2s_equal = isequal(E1, D2)
  
  
% % Differentiation of a hat function:
%   v = max(0,1-abs(x-pi)/2); clf
%   subplot(3,2,1), plot(x,v,'.-','markersize',13)
%   axis([0 2*pi -.5 1.5]), grid on, title('function')
%   subplot(3,2,2), plot(x,D*v,'.-','markersize',13)
%   axis([0 2*pi -1 1]), grid on, title('spectral derivative')
% 
% % Differentiation of exp(sin(x)):
%   v = exp(sin(x)); vprime = cos(x).*v;
%   subplot(3,2,3), plot(x,v,'.-','markersize',13)
%   axis([0 2*pi 0 3]), grid on
%   subplot(3,2,4), plot(x,D*v,'.-','markersize',13)
%   axis([0 2*pi -2 2]), grid on
%   error = norm(D*v-vprime,inf);
%   text(2.2,1.4,['max error = ' num2str(error)])