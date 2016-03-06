x = 0:.01:2*pi;
y = sin(x);
n = length(x);
k=0;
%for i=2:n
%    k=k+1;
yp=(sin(x(2:end))-sin(x(1:end-1)))./(x(2:end)-x(1:end-1));
%end
length(yp);
xp = x(2:end);
yp;
y_e =cos(xp);
length(xp);
length(y_e);
plot(xp,yp,xp,y_e)