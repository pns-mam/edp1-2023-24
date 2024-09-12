clear all, close all, clc
L = 5 ;            % longueur du domaine
nx = 100 ;         % nombre de mailles
dx = L/(nx+1) ;    % pas d'espace
cfl = 0.1 ;        % cfl
theta = 0.5;       % param du theta schema, =1/2 (Crank-Nicolson)
nu = 1;            % coefficient de diffusion
dt = dx^2/nu*cfl ; % dt = pas de temps
nt = 10 ;          % nombre de pas de temps effectues

x = linspace(0,L,nx+2); % maillage en espace
k = 2;
u0 = sin(k*pi/L*x); % solution initiale

alpha = (1-theta)*cfl; beta = theta*cfl;
A = diag((1+2*alpha)*ones(nx,1))+...
        diag(-alpha*ones(nx-1,1),1)+diag(-alpha*ones(nx-1,1),-1);
B = diag((1-2*beta)*ones(nx,1))+...
        diag(beta*ones(nx-1,1),1)+diag(beta*ones(nx-1,1),-1);
u = u0(2:end-1)' ;

for n = 1:nt
    u = A\(B*u) ;
    if mod(n,5) == 0
        plot(x,[0;u;0],'r-',x,u0,'b-'), grid on
        legend('schema implicite','donnee initiale')
        title('schema implicite');
        pause(0.1)
    end
end

% calcul de la solution exacte
uexacte = exp(-nu/L*(k*pi)^2*nt*dt)*sin(k*pi/L*x);

% comparaison au temps final
plot(x,[0;u;0],'r-',x,uexacte,'b-'), grid on
legend('theta schema','solution exacte')
title('theta schema');

% calcul de la norme de l'erreur 
norm([0;u;0]-uexacte','inf')
