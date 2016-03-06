% [x] = down_solve(L,b) takes a lower triangular matrix L and vector b
% and solves Lx=b

function [x] = down_solve(L,b)

[m,n] = size(L);

x(1) = b(1)/L(1,1);

for k=2:n
  x(k) = b(k);
  for j=1:k-1
    x(k) = x(k) - L(k,j)*x(j);
  end
  x(k) = x(k)/L(k,k);
end
x = x';