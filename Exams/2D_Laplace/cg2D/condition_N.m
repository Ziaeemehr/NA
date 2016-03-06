
nx = 4000;
ny = 40;
A = zeros(nx,ny);
for i = 1:nx
   for j = 1:ny
       if (i==j)
         A(i,j) = -4;
       end
       if ((i+1) == j) 
         A(i,j) = 1; 
       end
       if (i == (j+1)) 
         A(i,j) = 1;
       end
       if ((i+3) == j)
         A(i,j) = 1;
       end
       if (i == (j+3))
        A(i,j) = 1;       
       end
   end
end

condi=cond(A)