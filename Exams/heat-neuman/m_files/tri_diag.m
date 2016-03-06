% x = tri_diag(A,b) takes a tri-diagonal square matrix A and a vector b and 
% returns the solution of Ax=b.

function x = tri_diag(A,rhs)

[n,n] = size(A);

for i=1:n
  b(i) = A(i,i);
end
for i=1:n-1
  c(i) = A(i,i+1);
end
for i=2:n
  a(i) = A(i,i-1);
end

% solve for the entries of L and U so that LU = A:
beta(1) = b(1);
for j=2:n
  alpha(j) = a(j)/beta(j-1);
  beta(j) = b(j)-alpha(j)*c(j-1);
end

% solve Ly = b
y(1) = rhs(1);
for j=2:n
  y(j) = rhs(j)-alpha(j)*y(j-1);
end

% solve Ux = y
x(n) = y(n)/beta(n);
for j=1:n-1
  x(n-j) = (y(n-j)-c(n-j)*x(n-j+1))/beta(n-j);
end

x = x';
