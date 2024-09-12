function u = saddle(x,y)
u = (x-0.5)^2-(y-0.5)^2;
endfunction

function [A,b] = create_matrix(n,m)
dx = 1/(n-1);
dy = 1/(m-1);
dx2 = 1/dx^2;
dy2 = 1/dy^2;
x = 0:dx:1;
y = 0:dy:1;

A = spzeros(n*m,n*m);
b = zeros(n*m,1);
for i=1:n
 for j=1:m
   k = i+(j-1)*n;
   if (i==1 | i==n | j==1 | j==m) then
     b(k) = saddle(x(i),y(j));
     A(k,k) = 1.0;
   else
     b(k) = force;
     A(k,k)   = -2.0*(dx2+dy2);
     A(k,k+1) = dx2;
     A(k,k-1) = dx2;
     A(k,k+n) = dy2;
     A(k,k-n) = dy2;
   end
 end
end
endfunction


function [A,b] = create_matrix2(n,m)
dx = 1/(n-1);
dy = 1/(m-1);
dx2 = 1/dx^2;
dy2 = 1/dy^2;
x = 0:dx:1;
y = 0:dy:1;

d0 = -2*(dx2+dy2)*ones(n*m,1);
dmn = dy2*ones((m-1)*n,1);
dpn = dy2*ones((m-1)*n,1);
dm1 = dx2*ones(n*m-1,1);
dp1 = dx2*ones(n*m-1,1);

// conditions aux limites
d0(n:n:n*m) = 1;
d0(1:n:n*m) = 1;
d0(1:n) = 1;
d0(n*(m-1):n*m) = 1;

dmn(1:n:n*(m-1)) = 0;
dmn(n:n:n*(m-1)) = 0;
dmn(n*(m-2):n*(m-1)) = 0;

dpn(1:n:n*(m-1)) = 0;
dpn(n:n:n*(m-1)) = 0;
dpn(1:n) = 0;

dm1(n-1:n:n*m-1) = 0;
dm1(n:n:n*m-1) = 0;
dm1(1:n-1) = 0;
dm1(n*(m-1):n*m-1) = 0;

dp1(n+1:n:n*m-1) = 0;
dp1(n:n:n*m-1) = 0;
dp1(1:n-1) = 0;
dp1(n*(m-1):n*m-1) = 0;

A = diag(sparse(d0))+diag(sparse(dm1),-1)+diag(sparse(dp1),1)+...
diag(sparse(dmn),-n)+diag(sparse(dpn),n);

b = force*ones(n*m,1);
// conditions aux limites

b(1:n)           = feval(x,y(1),saddle);
b(n*(m-1)+1:n*m) = feval(x,y(m),saddle); 
b(1:n:n*(m-1)+1) = feval(x(1),y,saddle)';
b(n:n:n*m)       = feval(x(n),y,saddle)'; 

endfunction


function spy(A)
[i,j] = find(A~=0);
[N,M] = size(A);
xsetech([0,0,1,1],[1,0,M+1,N])
xrects([j;N-i+1;ones(i);ones(i)],ones(i))
xrect(1,N,M,N)
endfunction

function U = solve_equation(A,b)
  U = A\b;
endfunction

function run_problem(n,m)

dx = 1/(n-1);
dy = 1/(m-1);
x = 0:dx:1;
y = 0:dy:1; 

printf('Create Matrix:')
timer()
[A,b] = create_matrix2(n,m);
printf('time = %g\n',timer())

printf('Solve Matrix problem:')
timer()
U = solve_equation(A,b);
printf('time = %g\n',timer())

u = matrix(U,n,m);

xbasc();
xselect();
plot3d1(x,y,u)
endfunction


