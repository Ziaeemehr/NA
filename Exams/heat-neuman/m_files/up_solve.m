% [x] = up_solve(U,b) takes an upper triangular matrix U and vector b
% and solves Ux=b

function [x] = up_solve(U,b)

[m,n] = size(U);

x(n) = b(n)/U(n,n);

for k=1:n-1
  x(n-k) = b(n-k);
  for j=n-k+1:n
    x(n-k) = x(n-k) - U(n-k,j)*x(j);
  end
  x(n-k) = x(n-k)/U(n-k,n-k);
end
x = x';
