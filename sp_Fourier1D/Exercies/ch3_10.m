% p6.m - variable coefficient wave equation

No_ = 10;

T = 13; %- This should be checked by analytic solution of cirularity
Nx = [4, 5, 6, 7, 8];%, 9, 10, 11, 12, 13];  %-- The power of 2 : we can expenad this independently.
ns = size(Nx, 2);
Tdx = zeros(1, ns);   %   array of maximum difference between v0 and vT


for in = 1:ns

% Grid, variable coefficient, and initial data:
  N = 2^Nx(in); 

  h = 2*pi/N; 
  x = h*(1:N); t = 0; dt = h/4;
  c = .2 + sin(x-1).^2;
  v = exp(-100*(x-1).^2); vold = exp(-100*(x-.2*dt-1).^2);

% Time-stepping by leap frog formula:
  tplot = .15; clf, drawnow, set(gcf,'renderer','zbuffer')
  plotgap = round(tplot/dt); dt = tplot/plotgap;
  nplots = round(T/tplot);
  data = [v; zeros(nplots,N)]; tdata = t;
  for i = 1:nplots
    for n = 1:plotgap
      t = t+dt;
      v_hat = fft(v);
      w_hat = 1i*[0:N/2-1 0 -N/2+1:-1] .* v_hat;
      w = real(ifft(w_hat)); 
      vnew = vold - 2*dt*c.*w; vold = v; v = vnew;
    end
    data(i+1,:) = v; tdata = [tdata; t];
    
    if i == nplots
        v0 = data(1, :)
        vT = data(nplots+1, :)
    end
  end

  maxdif = max(abs(v0-vT));
  
  Tdx(in) = maxdif;
  
end

plot(Nx, Tdx, 'r');
xlabel('log(2)N');
ylabel('maximum differences');
  for i = 1:ns
      text(Nx(i), Tdx(i), num2str(Tdx(1, i)));
  end

%   waterfall(x,tdata,data), view(10,70), colormap(1e-6*[1 1 1]);
%   axis([0 2*pi 0 T 0 5]), ylabel t, zlabel u, grid off
