// Matrice du laplacien en differences finies
function A = laplaceD(n)
A = n^2*(diag(2*ones(n-1,1))+diag(-ones(n-2,1),-1)+diag(-ones(n-2,1),1));
endfunction

// Inialisation du second membre
function b = initrhs(n,func) 
dx = 1/n;
x = linspace(dx,1-dx,n-1);
b = evstr(func)';
endfunction

// Methode de Jacobi
function [x,iter] = Jacobi(A,b,tol,iterMax,x)

[sortie, entree] = argn(0);
if entree == 4 then, x = zeros(b); end;
if entree == 3 then, x = zeros(b); iterMax = 200; end;
if entree == 2 then, x = zeros(b); iterMax = 200; tol = 1.0e-4; end;

D = diag(1 ./diag(A));

iter = 0; r = b-A*x;
while (norm(r) > tol) & (iter < iterMax) do
  x = x + D*r;
  r = b-A*x;
  iter = iter + 1;
end;

endfunction

// Methode de relaxation
function [x,iter] = Relax(A,b,omega,tol,iterMax,x)

[sortie, entree] = argn(0);
if entree == 5 then, x = zeros(b); end;
if entree == 4 then, x = zeros(b); iterMax = 200; end;
if entree == 3 then, x = zeros(b); iterMax = 200; tol = 1.0e-4; end;
if entree == 2 then, x = zeros(b); iterMax = 200; tol = 1.0e-4; omega = 1; end;

M = diag((1-omega)/omega*diag(A)) + tril(A);

iter = 0; r = b-A*x;
while (norm(r) > tol) & (iter < iterMax) do
  x = x + M\r;
  r = b-A*x;
  iter = iter + 1;
end;

endfunction

// Algorithme du gradient a pas constant
function [x,iter] = Gradient(A,b,alpha,tol,iterMax,x)

[sortie, entree] = argn(0);
if entree == 5 then, x = zeros(b); end;
if entree == 4 then, x = zeros(b); iterMax = 2000; end;
if entree == 3 then, x = zeros(b); iterMax = 2000; tol = 1.0e-4; end;
if entree == 2 then, x = zeros(b); iterMax = 2000; tol = 1.0e-4; alpha = 1.0e-4; end;

iter = 0;r = b-A*x;
while (norm(r) > tol) & (iter < iterMax) do
  x = x+alpha*r;
  r = b-A*x;
  iter = iter +1;
end;
endfunction

// Algorithme du gradient a pas variable
function [x,iter] = GradientV(A,b,tol,iterMax,x)

[sortie, entree] = argn(0);
if entree == 5 then, x = zeros(b); end;
if entree == 4 then, x = zeros(b); iterMax = 2000; end;
if entree == 3 then, x = zeros(b); iterMax = 2000; tol = 1.0e-4; end;
if entree == 2 then, x = zeros(b); iterMax = 2000; tol = 1.0e-4; alpha = 1.0e-4; end;

iter = 0;r = b-A*x;
alpha = norm(r)/(r'*A*r);
while (norm(r) > tol) & (iter < iterMax) do
  x = x+alpha*r;
  r = b-A*x;
  iter = iter +1;
  alpha = norm(r)/(r'*A*r);
end;
endfunction

// TEST POUR LA METHODE DE JACOBI

n  = 20; A  = laplaceD(n);
dx = 1/n; x = linspace(dx,1-dx,n-1);
b  = initrhs(n,'x.*sin(x)');
sol = A\b;
tol  = 1.0e-2; x=Jacobi(A,b,tol,1000); norm(x-sol), norm(inv(A))*tol 
tol  = 1.0e-3; x=Jacobi(A,b,tol,1000); norm(x-sol), norm(inv(A))*tol 
tol  = 1.0e-4; x=Jacobi(A,b,tol,1000); norm(x-sol), norm(inv(A))*tol 

// TEST POUR LA METHODE DE RELAXATION

n  = 20; A  = laplaceD(n);
dx = 1/n; x = linspace(dx,1-dx,n-1);
b  = initrhs(n,'sin(x)');
sol = A\b;

step = 0.1;
omega = 0.1:0.1:2;
for i = 1:length(omega)
  [x,iter] = Relax(A,b,omega(i),1.0e-6,1000);
  it(i) = iter;
end;

plot(omega,it,'bo-');


// TEST POUR LA METHODE DU GRADIENT

n  = 10; A  = laplaceD(n);
dx = 1/n; x = linspace(dx,1-dx,n-1);
b  = initrhs(n,'x.*sin(x)');  
[x,iter] = Gradient(A,b,1.0e-4,1.0e-4,10000);
norm(x-A\b)

alpha = 40*1.0e-4:1.0e-5:60*1.0e-4;
for i = 1:length(alpha)
  [x,iter] = Gradient(A,b,alpha(i),1.0e-10,2000);
  nit(i) = iter;
end;
clf()
plot(alpha,nit)

// TEST POUR LA METHODE DU GRADIENT A PAS VARIABLE

n  = 5; A  = laplaceD(n);
dx = 1/n; x = linspace(dx,1-dx,n-1);
b  = initrhs(n,'x.*sin(x)');  
lam = spec(A);alphaopt = 2/(min(lam)+max(lam));
[x1,iter1] = Gradient(A,b,alphaopt,1.0e-4,10000);
[x1,iter2] = GradientV(A,b,1.0e-4,10000);
// on doit observer que iter1 >> iter2
