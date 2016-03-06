% [u,err,x,t] = heat3(t_0,t_f,M,N)
%
% this solves the heat equation u_t = u_xx with initial data u_0 =
% cos(x) with periodic boundary conditions using finite-differences in
% space and implicit time-stepping.  t_0 is the initial time, t_f is the
% final time, N is the number of mesh-points, and M is the number of
% time steps.  err is the error. 

function [u,err,x,t] = heat2(t_0,t_f,M,N)

% define the mesh in space
dx = 2*pi/(N-1);
x = 0:dx:2*pi-dx;
x = x';

% define the mesh in time
dt = (t_f-t_0)/M;
t = t_0:dt:t_f;

% define the ratio r
r = dt/dx^2 

for i=1:N-1
  u(i,1) = cos(x(i));
end 
err(:,1) = u(:,1) - exp(t_0-t(1))*cos(x);

% define the matrix A which has to be inverted at every time-step.
%   u_new(1) - u_old(1) = dt/h^2 (u_new(2)-2*u_new(1)+u_new(N-1)
%   u_new(i) - u_old(i) = dt/h^2 (u_new(i+1)-2*u_new(i)+u_new(i-1)
%   u_new(N) - u_old(N) = dt/h^2 (u_new(2)-2*u_new(N)+u_new(N-1)

A = zeros(N-1,N-1);
A(1,1) = 1+2*r;
A(1,2) = -r;
A(1,N-1) = -r;
for i=2:N-2
  A(i,i-1) = -r;
  A(i,i) = 1+2*r;
  A(i,i+1) = -r;
end
A(N-1,N-2) = -r;
A(N-1,N-1) = 1+2*r;
A(N-1,1) = -r;
% note: I'm not keeping around the Nth meshpoint since x(N)=x(1).  This has the effect
% that A is tri-diagonal matrix except with scruff in the top right and bottom left
% corner.  I'll then use a subroutine periodic_tridiag to invert this matrix
%
% as we saw, using the LU approach was impossibly slow, so I really do want the O(n)
% option, even if it meant that I had to write an extra subroutine.

for j=1:M
  [u(:,j+1)] = periodic_tridiag(A,u(:,j));
  err(:,j+1) = u(:,j+1) - exp(t_0-t(j+1))*cos(x);
end

u(N,:) = u(1,:);
err(N,:) = err(1,:);
x(N) = 2*pi;
