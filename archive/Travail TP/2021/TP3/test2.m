clear all, close all, clc
% Schema de Lax-Wendroff
V = 0.1;     % vitesse d’advection
Tmax = 10;   % temps maximum
L = 5;       % longueur du domaine

N = 1001;      % nb de points de discretisation en espace
dx = L/(N-1); % pas d’espace
s = 0.9;      % nombre de Courant
dt = s*dx/V;  % pas de temps

x = linspace(0,L,N)';
condinit = @(x) (x>1.0) & (x<1.5);
%condinit = @(x) sin(8*pi*x/L);
u = double(condinit(x));

tps = dt:dt:Tmax;
for t = tps
    uold = u;
    uexact = condinit(x-V*t);
    u(N) = uexact(N);
    u(2:N-1) = uold(2:N-1)-s*(uold(3:N)-uold(1:N-2))/2 +...
        s^2*(uold(3:N)-2*uold(2:N-1)+uold(1:N-2))/2 ;
    u(1) = u(N);
end

plot(x,u,'rx-',x,uexact,'bo-'), grid on
title('Comparaison entre la solution exacte et approchee')
legend('Sol. approchee','Sol. exacte')
xlabel('x')
ylabel('Solution')