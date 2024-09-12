clear all, close all, clc
L = 5 ;            % longueur du domaine
nx = 20 ;          % nombre de mailles
dx = L/(nx+1) ;    % pas d'espace
cfl = 0.01 ;       % cfl
nu = 1;            % coefficient de diffusion
dt = dx^2/nu*cfl ; % dt = pas de temps
nt = 20 ;          % nombre de pas de temps effectues

x = linspace(0,L,nx+2); % maillage en espace
k = 5;
u0 = sin(k*pi/L*x); % solution initiale

A = diag((1+2*cfl)*ones(nx,1))+diag(-cfl*ones(nx-1,1),1)+diag(-cfl*ones(nx-1,1),-1);
u = u0(2:end-1)' ;

for n = 1:nt
    u = A\u ;
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
legend('schema implicite','solution exacte')
title('schema implicite');

% calcul de la norme de l'erreur 
norm([0;u;0]-uexacte','inf')