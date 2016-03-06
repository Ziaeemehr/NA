% [u,err,x,t] = heat2(t_0,t_f,M,N)
%
% this solves the heat equation u_t = u_xx with initial data u_0 =
% cos(x) with periodic boundary conditions using finite-differences in
% space and implicit time-stepping.  t_0 is the initial time, t_f is the
% final time, N is the number of mesh-points, and M is the number of
% time steps.  err is the error. 

function [u,err,x,t] = heat2(t_0,t_f,M,N)

% define the mesh in space
dx = 2*pi/(N-1);
x = 0:dx:2*pi;
x = x';

% define the mesh in time
dt = (t_f-t_0)/M;
t = t_0:dt:t_f;

% define the ratio r
r = dt/dx^2 

for i=1:N
  u(i,1) = cos(x(i));
end 
err(:,1) = u(:,1) - exp(t_0-t(1))*cos(x);

% define the matrix A which has to be inverted at every time-step.
%   u_new(1) - u_old(1) = dt/h^2 (u_new(2)-2*u_new(1)+u_new(N-1)
%   u_new(i) - u_old(i) = dt/h^2 (u_new(i+1)-2*u_new(i)+u_new(i-1)
%   u_new(N) - u_old(N) = dt/h^2 (u_new(2)-2*u_new(N)+u_new(N-1)
A = zeros(N,N);
A(1,1) = 1+2*r;
A(1,2) = -r;
A(1,N-1) = -r;
for i=2:N-1
  A(i,i-1) = -r;
  A(i,i) = 1+2*r;
  A(i,i+1) = -r;
end
A(N,2) = -r;
A(N,N-1) = -r;
A(N,N) = 1+2*r;

% lazy soln: do an LU decomposition on A and then solve A u_new = u_old.
% fast soln: figure out how to do the corners so that I can then use a tri-diagonal
% solver.
%
[L,UU] = my_LU(A);

% note: matlab can tell the difference between u and U, so I could have U as the  
% matrix and u as the solution.  But not all languages will do this, sometimes you
% have to tell them explicitly to care about case.  The safest thing is to not have the
% habit of depending on cases.

% then we'll just do
%[y] = down_solve(L,u_old);
%[u_new] = up_solve(UU,y);

for j=1:M
  [y] = down_solve(L,u(:,j)); 
  [u(:,j+1)] = up_solve(UU,y);
  err(:,j+1) = u(:,j+1) - exp(t_0-t(j+1))*cos(x);
end
% note: in the above, we really took advantage of the fact that matlab
% is a matrix-driven language.  Otherwise we would have had to do an
% extra for loop on the i-index.

