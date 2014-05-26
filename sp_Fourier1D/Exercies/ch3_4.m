% ch3_4.m - Question that shows the difference in error between Differentation matrix method and method using FFT
%           for the function v = exp(sin(x))

% For complex v, delete "real" commands.
% Set up grid and differentiation matrix:

No_ = 4

vmin = 24;
vmax = 100;
Nv = vmax-vmin;
Nu = Nv/2 + 1;
w0 = [vmin : 2 : vmax]';
wd = zeros(Nu, 1);
wf = zeros(Nu, 1);

i = 1;
for iN = vmin : 2 : vmax

  h = 2*pi/iN; 
  x = h*(1:iN)';

% Differentiation of exp(sin(x)):
  v = exp(sin(x)); 
  vprime = cos(x) .* v;
 
% 1. Spectral Differentiation using matrix
  column = [0 .5*(-1).^(1:iN-1).*cot((1:iN-1)*h/2)]';
  D = toeplitz(column,column([1 iN:-1:2]));
  error_D = norm(D*v-vprime,inf);
  
% 2. Spectral Differentiation using FFT
%    For complex v, delete "real" commands.
  v_hat = fft(v);
  w_hat = 1i .* [0:iN/2-1 0 -iN/2+1:-1]' .* v_hat;
  w = real(ifft(w_hat));
  error_FFT = norm(w-vprime,inf);
  
% 3. Store errors to wd, wf for display
  wd(i) = error_D;
  wf(i) = error_FFT;
  i = i  + 1;

end

%Subplot 1
  subplot(1,2,1);
  plot(w0,wd,'b', w0,wf,'r');
  legend('Er-Diff. Matrix','Er-FFT',0);
  %axis(0, 1e-14, 0, 100);
  grid on   
  xlabel('N');
  title('Errors');
  
%Subplot 2
  subplot(1,2,2);
  plot(w0,(wd-wf)*100/wd);
  legend('Given by (wd-wf)*100/wd',0);
  grid on;
  xlabel('N');
  ylabel('%');
  title('Error difference (%)');


