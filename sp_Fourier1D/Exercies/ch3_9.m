% ch3_9.m - variable coefficient wave equation

No_ = 9;
Nx = [128, 160, 192, 224, 256];%, 384, 512];  -- we can expenad this independently.
ns = size(Nx, 2);
Tmx = zeros(1, ns);   %   array of elapsed time for Finite Difference Matrix
Tfx = zeros(1, ns);   %   array of elapsed time for FFT

for in = 1:ns

% Grid, variable coefficient, and initial data:
  N = Nx(in); 

% 1. Finite difference spectral derivative
% Grid, variable coefficient, and initial data:
  h = 2*pi/N; 
  x = h*(1:N); 
  t = 0; 
  dt = h/4;
  c = .2 + sin(x-1).^2;
  v = exp(-100*(x-1).^2);
  vold = exp(-100*(x-.2*dt-1).^2);

% Time-stepping by leap frog formula:
  tmax = 8; tplot = .15; clf, drawnow, set(gcf,'renderer','zbuffer')
  plotgap = round(tplot/dt); dt = tplot/plotgap;
  nplots = round(tmax/tplot);
  
  datad = [v; zeros(nplots,N)]; 
  tdatad = t;
  
  tic;
  
  for i = 1:nplots
    for n = 1:plotgap
      t = t+dt;
     
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      x1 = [0 x h*(N+1)];
      v1 = interp1(x, v, x1, 'spline');
      vf = v1(3:N+2); vb = v1(1:N);
      w = (vf - vb)./(2*h);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      vnew = vold - 2*dt*c.*w; vold = v; v = vnew;
    end
    datad(i+1,:) = v; tdatad = [tdatad; t];
  end

  ti = toc;
  Tmx(in) = ti;


% 2. FFT spectral derivative
% Grid, variable coefficient, and initial data:
  h = 2*pi/N; 
  x = h*(1:N); 
  t = 0; 
  dt = h/4;
  c = .2 + sin(x-1).^2;
  v = exp(-100*(x-1).^2);
  vold = exp(-100*(x-.2*dt-1).^2);

% Time-stepping by leap frog formula:
  tmax = 8; tplot = .15; clf, drawnow, set(gcf,'renderer','zbuffer')
  plotgap = round(tplot/dt); dt = tplot/plotgap;
  nplots = round(tmax/tplot);

  dataf = [v; zeros(nplots,N)]; 
  tdataf = t;  
  
  tic;
  
  for i = 1:nplots
    for n = 1:plotgap
      t = t+dt;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      v_hat = fft(v);
      w_hat = 1i*[0:N/2-1 0 -N/2+1:-1] .* v_hat;
      w = real(ifft(w_hat));       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      vnew = vold - 2*dt*c.*w; vold = v; v = vnew;
    end
    dataf(i+1,:) = v; tdataf = [tdataf; t];
  end
  
  tj = toc;
  Tfx(in) = tj;
  
end


figure(1);

  plot(Nx, Tmx, 'r', Nx, Tfx, 'g');
  legend('T: Finite Difference','T: FFT', 0);
  grid on   
  xlabel('N');
  title('Elapsed Time Comparison between Finite Difference, FFT for Sp.Diff.');
  for i = 1:ns
      text(Nx(i), Tmx(i), num2str(Tmx(1, i)));
      text(Nx(i), Tfx(i), num2str(Tfx(1, i)));
  end

  
figure(2);

subplot(1, 2, 1);

  waterfall(x, tdatad, datad), view(10,70), colormap(1e-6*[1 1 1]);
  axis([0 2*pi 0 tmax 0 5]);
  ylabel t, zlabel u, grid off
  title('Finite Difference');

subplot(1, 2, 2);

  waterfall(x, tdataf, dataf), view(10,70), colormap(1e-6*[1 1 1]);
  axis([0 2*pi 0 tmax 0 5]);
  ylabel t, zlabel u, grid off
  title('FFT');
