% [u,err,x,t,kk,amp] = heat4(t_0,t_f,M,N)
%
% this solves the periodic heat equation u_t = u_xx with initial data u_0 =
% cos(x) with periodic boundary conditions using spectral derivatives in
% space and explicit time-stepping.  t_0 is the initial time, t_f is the
% final time, N is the number of mesh-points, and M is the number of
% time steps.  err is the error. it also returns the power spectrum of the
% solution, where kk is the collection of wave numbers and amp = log(f_k^2)
% is the magnitude of the kth Fourier mode.
%
% N SHOULD BE A POWER OF 2!  
%
% (It will work no matter what, but will go much faster if N is a power
% of 2.)

function [u,err,x,t,kk,amp] = heat4(t_0,t_f,M,N)

% define the mesh in space
dx = 2*pi/N;
x = 0:dx:2*pi-dx;
x = x';
kk = 0:1:(N/2-1);

% define the mesh in time
dt = (t_f-t_0)/M;
t = t_0:dt:t_f;

% define the ratio r
r = dt/dx^2 
% define the maximum amplification
r = 1 - dt*(N/2-1)^2

for i=1:N
  u(i,1) = cos(x(i));
end 
err(:,1) = u(:,1) - exp(t_0-t(1))*cos(x);
v = fft(u(:,1));
for k=0:(N/2-1)
  if k == 0 
    amp(k+1,1) = log10(sqrt(v(1+k)^2)/N + 10e-16);
  else
    amp(k+1,1) = log10(sqrt(real(v(1+k)*v(N-k+1)))/N + 1e-16);    
  end
end

% take the fourier transform of u_old.  At each mode, multiply by 1 -
% dt*k^2 where k is the mode-number.  Then take the inverse fourier
% transform

for j=1:M
  v = fft(u(:,j));
  for k = 1:(N/2-1)
% k is the wave-number.  Note that the zero mode v(1) is untouched since
% it is the mean.
    v(1+k) = (1-dt*k^2)*v(1+k);
    v(N-k+1) = (1-dt*(-k)^2)*v(N-k+1);
  end
  for k=0:(N/2-1)
    if k == 0 
      amp(k+1,j+1) = log10(v(1+k)/N + 10e-16);
    else
      amp(k+1,j+1) = log10(sqrt(real(v(1+k)*v(N-k+1)))/N + 10e-16);    
    end
  end
% take the inverse transform.  I've cut off the imaginary parts since
% they're going to be at the level of round-off anyway.  (the true
% solution started out real-valued and will remain real.  If I wanted to
% be super-cautious, I'd keep the imaginary parts of the computed
% solution around to make sure they aren't growing on me.)
  u(:,j+1) = real(ifft(v)); 
  err(:,j+1) = u(:,j+1) - exp(t_0-t(j+1))*cos(x);
end

