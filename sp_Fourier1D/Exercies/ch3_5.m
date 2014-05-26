% ch3_5.m - Question that capturing calculation time of FFT

No_ = 5

vmin = 1;
vmax = 15;
repet = 10;
Nu = vmax-vmin+1;
w0 = [vmin : 1 : vmax]';
wt = zeros(Nu, 1);

for j = vmin : 2 : vmax

  iN = 2^j;
  h = 2*pi/iN; 
  x = h*(1:iN)';

% Differentiation of exp(sin(x)):
  v = exp(sin(x)); 
   
% Spectral Differentiation using FFT
% For complex v, delete "real" commands.

tic;

  for r=1:repet
    v_hat = fft(v);
    w_hat = 1i .* [0:iN/2-1 0 -iN/2+1:-1]' .* v_hat;
    w = real(ifft(w_hat));
  end

t = toc;
  wt(j) = t;
  
end

%Subplot 1
  subplot(1,2,1);
  plot(w0, wt, 'r');
  legend('Time elapsed for FFT',0);
  grid on   
  xlabel('N:The Power of 2');
  title('FFT of 100 times computation');
  
  
  
%==========================================================================

vmin = 500;
vmax = 520;
% vmin = 128;
% vmax = 1024;

repet = 10;
Nv = vmax-vmin+1;
%Nv = vmax-vmin;
Nu = Nv;
%Nu = Nv/2 + 1;

w0 = [vmin : 1 : vmax]';
% w0 = [vmin : 2 : vmax]';

wt = zeros(Nu, 1);

i = 1;
for iN = vmin : 2 : vmax

  h = 2*pi/iN; 
  x = h*(1:iN)';

% Differentiation of exp(sin(x)):
  v = exp(sin(x)); 
   
% Spectral Differentiation using FFT
% For complex v, delete "real" commands.

tic;

  for r=1:repet
    v_hat = fft(v);
    w_hat = 1i .* [0:iN/2-1 0 -iN/2+1:-1]' .* v_hat;
    w = real(ifft(w_hat));
  end

t = toc;
  wt(i) = t;
  i = i + 1;
  
end

%-------------------------------------------------------------------------------------
% % Prime element testing
%   nP = 0;
%   for iP = vmin : 1 : vmax
%     if isprime(iP) == 1
%       nP = nP + 1;
%     end
%   end
%   
%   % Matrices of primes wp[] and elapsed times wpt[]
%   i = 1;
% 
%   for iP = vmin : 1 : vmax
%     if isprime(iP) == 1
%       wp(i) = iP;
%       
%         h = 2*pi/iP; 
%         x = h*(1:iP)';
% 
%         % Differentiation of exp(sin(x)):
%         v = exp(sin(x)); 
%    
%       tic;
% 
%         for r=1:repet
%           v_hat = fft(v);
%                a = size(v_hat)
%                b =size([0:iP/2-1 0 -iP/2+1:-1]')
%         w_hat = 1i .* [0:iP/2-1 0 -iP/2+1:-1]' .* v_hat;
%           b = size(w_hat)
%           
%           w = real(ifft(w_hat));
%         end
% 
%       t = toc
%       wpt(i) = t;
% 
%       i = i + 1;
%     end
%   end
%  
% 
  %Subplot 2 
  subplot(1,2,2);
%  plot(w0, wt, 'g', wp, wpt, 'r');
  plot(w0, wt, 'g');
  legend('Time elapsed for FFT',0);
  grid on   
  xlabel('N');
  title('FFT computation');