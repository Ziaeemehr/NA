% Iterative Solutions of linear euations:(1) Jacobi Method
% Linear system: A x = b
% Coefficient matrix A, right-hand side vector b
A=[7 3 -1 2; 3 8 1 -4; -1 1 4 -1; 2 -4 -1 6];
b= [-1;0;-3;1];
% Set initial value of x to zero column vector 
x0=zeros(1,4);
% Set Maximum iteration number k_max
k_max=1000;
% Set the convergence control parameter erp
erp=0.0001;
% Show the q matrix
%q=diag(diag(A))
% loop for iterations
for k=1:k_max
   for i=1:4
      s=0.0;
      for j=1:4
         if j==i 
             continue
         else    
             s=s+A(i,j)*x0(j);
         end
      end
      x1(i)=(b(i)-s)/A(i,i);
   end
   if norm(x1-x0)<erp
      break
   else
      x0=x1;   
   end
end
% show the final solution
x=x1
% show the total iteration number
n_iteration=k

