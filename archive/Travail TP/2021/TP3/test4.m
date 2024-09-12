clear all, close all, clc
% Schema de Lax-Wendroff implicite
V = 0.1;     % vitesse d’advection
Tmax = 10;   % temps maximum
L = 5;       % longueur du domaine

N = 1001;      % nb de points de discretisation en espace
dx = L/(N-1); % pas d’espace
s = 0.9;      % nombre de Courant
dt = s*dx/V;  % pas de temps

x = linspace(0,L,N)';
cinit = 1;    % 1 = creneau, 2 = sinus

x = linspace(0,L,N)';
if cinit == 1
    condinit = @(x) (x>1.0) & (x<1.5);
else
    condinit = @(x) sin(8*pi*x/L);
end
u = double(condinit(x));

A = diag((1+s^2)*ones(N,1))-diag(s^2/2*ones(N-1,1),-1)-...
    diag(s^2/2*ones(N-1,1),1); 
A(1,N) = -s^2/2; % on modifie la premiere ligne et la derniere
A(N,1) = -s^2/2; % afin de prendre en compte la periodicite

tps = dt:dt:Tmax;
w = zeros(size(u));

[L,U,P]= lu(A); % on calcule la factorisation LU de A
tic
for t = tps
    uold = u;
    uexact = condinit(x-V*t);
    % second membre dans le schema implicite
    w(N) = uexact(N); 
    w(2:N-1) = uold(2:N-1)-s/2*(uold(3:N)-uold(1:N-2));
    w(1) = w(N);
    %resolution du systeme lineaire
    u = U\(L\(P*w));
    plot(x,u,'rx-',x,uexact,'bo-'), grid on
    title('Comparaison entre la solution exacte et approchee')
    legend('Sol. approchee','Sol. exacte')
    xlabel('x')
    ylabel('Solution')
    %pause(0.1)
end
toc
