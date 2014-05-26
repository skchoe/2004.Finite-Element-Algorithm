
%--------------------------------------------------------------------------------------
function [w, etime] = spectraldiffbyfft(v, N, h, bt);
  
  if bt ~= 0
      tic;
  end

  
  v_hat = fft(v);
  w_hat = 1i*[0:N/2-1 0 -N/2+1:-1] .* v_hat;
  w = real(ifft(w_hat));
    
  if bt ~= 0
    etime = toc;
  else 
    etime = 0;
  end

  
return
  