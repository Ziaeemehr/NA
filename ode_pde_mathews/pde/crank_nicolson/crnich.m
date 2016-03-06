function U = crnich(c1,c2,a,b,c,n,m)

h = a/(n-1);
k = b/(m-1);
r = c*c *k/(h*h);
s1 = 2 + 2/r;
s2 = 2/r - 2;
U = zeros(n,m);
% B.C
U(1,1:m) = c1;
U(n,1:m) = c2;
xx = h:h:(n-2)*h;
%U(2:n-1,1) = feval(f,h:h:(n-2)*h);
U(2:n-1,1) = sin(pi*xx) + sin(3*pi*xx);
Vd(1,1:n) = s1*ones(1,n);
Vd(1)=1;
Vd(n)=1;
Va = -ones(1,n-1);
Va(n-1)=0;
Vc =-ones(1,n-1);
Vc(1) = 0;
Vb(1) = c1;
Vb(n) = c2;
for i=2:n
    for j = 2:n
        if (i==j)     
            A(i,j) = s1;           
        elseif (i==(j+1)) 
            A(i,j) = -1;
        elseif ((i+1)==j) 
            A(i,j) = -1;
        end
    end
end
A(n,n)=1;      A(1,2) = 0;
A(n,n-1) = 0;  A(1,1) = 1;

for j = 2:m
    for i = 2:n-1
        Vb(i) = U(i-1,j-1) + U(i+1,j-1) + s2 * U(i,j-1);
    end
     X = A\Vb';
    U(1:n,j) = X';
end
U = U';