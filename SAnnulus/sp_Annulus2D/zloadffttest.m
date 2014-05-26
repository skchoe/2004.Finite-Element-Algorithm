zloadffttest = 0
% 
% N = 14;
% h = 2*pi/N;
% x = h * (0:1:N-1)';
% x = [x']';
% 
% x1 = cos(x);
% %x1 = x.*x+1;
% t = fft(x1)
% 
% 

N = 14;
h = 2*pi/N;
u = h * (0:1:N-1)'

w = fft(u)
g = ifft(w)