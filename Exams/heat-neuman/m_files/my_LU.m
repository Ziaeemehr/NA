% [L,U] = LU(A) LU takes a square matrix A and returns upper and
% lower triangular matrices such that LU = A.  NOTE: it writes over A!
%
function [L,U] = LU(A);

[m,n] = size(A);

% note that I'm not initializing the entire matrix L and U since it
% turns out that Matlab will assume things are zero automatically. 
for i=1:n
  L(i,i) = 1;
end

for k=1:n-1
  for i=k+1:n
    mult = A(i,k)/A(k,k);
    for j=k+1:n
      A(i,j) = A(i,j) - A(k,j)*mult;
    end
    L(i,k) = mult;
  end
end

for i=1:n
  for j=i:n
    U(i,j) = A(i,j);
  end
end
