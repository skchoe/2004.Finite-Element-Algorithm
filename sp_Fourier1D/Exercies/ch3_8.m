% ch3_8.m - Looking at two Time stepping method: leap frog/matlab built-in 'expm'
No_ = 8

Nx = [128, 160, 192, 224, 256];%, 384, 512];  -- we can expenad this independently.
ns = size(Nx, 2);
Tmx = zeros(1, ns);   %   array of elapsed time for spectral differentiation matrix
Tfx = zeros(1, ns);   %   array of elapsed time for FFT

for in = 1:ns

% Grid, variable coefficient, and initial data:
  N = Nx(in); 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For FFT differentiation
  h = 2*pi/N; x = h*(1:N); t = 0; dt = h/4;
  c = .2 + sin(x-1).^2;
  v = exp(-100*(x-1).^2); 
  v_old = exp(-100*(x-.2*dt-1).^2); % init for v(-1)

  % Time-stepping by leap frog formula:
  tmax = 8; tplot = .15; clf, drawnow, set(gcf,'renderer','zbuffer')
  plotgap = round(tplot/dt);
  dt = tplot/plotgap;
  nplots = round(tmax/tplot);
  datam = [v; zeros(nplots,N)]; 
  tdatam = t;

  % Spectral differentiation matrix
  tic;
  
  for i = 1:nplots
    for n = 1:plotgap
      t = t + dt;
   
      [w, etime] = spectraldiffbymatrix(v, N, h, 0);
      
      v_new = v_old - 2*dt*c.*w; 
      v_old = v; 
      v = v_new;
    end
    
    datam(i+1,:) = v; 
    tdatam = [tdatam; t];
  end
  
  ti = toc;
  Tmx(in) = ti;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For FFT differentiation
  h = 2*pi/N; x = h*(1:N); t = 0; dt = h/4;
  c = .2 + sin(x-1).^2;
  v = exp(-100*(x-1).^2); 
  v_old = exp(-100*(x-.2*dt-1).^2); % init for v(-1)

  % Time-stepping by leap frog formula:
  tmax = 8; tplot = .15; clf, drawnow, set(gcf,'renderer','zbuffer')
  plotgap = round(tplot/dt);
  dt = tplot/plotgap;
  nplots = round(tmax/tplot);
  dataf = [v; zeros(nplots,N)]; 
  tdataf = t;

  tic;
  % FFT
  for i = 1:nplots
    for n = 1:plotgap
      t = t + dt;
   
      [w, etime] = spectraldiffbyfft(v, N, h, 0);
      
      v_new = v_old - 2*dt*c.*w; 
      v_old = v; 
      v = v_new;
    end
    
    dataf(i+1,:) = v; 
    tdataf = [tdataf; t];
  end
  
  tj = toc;
  Tfx(in) = tj;
  
end %End of all loop

figure(1);
  plot(Nx, Tmx, 'r', Nx, Tfx, 'g');
  legend('T: SP. Diff Matrix','T: FFT', 0);
  grid on   
  xlabel('N');
  title('Elapsed Time Comparison between Matrix, FFT for Sp.Diff.');
  for i = 1:ns
      text(Nx(i), Tmx(i), num2str(Tmx(1, i)));
      text(Nx(i), Tfx(i), num2str(Tfx(1, i)));
  end

  
  
figure(2);

subplot(1, 2, 1);

  waterfall(x, tdatam, datam), view(10,70), colormap(1e-6*[1 1 1]);
  axis([0 2*pi 0 tmax 0 5]);
  ylabel t, zlabel u, grid off

subplot(1, 2, 2);

  waterfall(x, tdataf, dataf), view(10,70), colormap(1e-6*[1 1 1]);
  axis([0 2*pi 0 tmax 0 5]);
  ylabel t, zlabel u, grid off
