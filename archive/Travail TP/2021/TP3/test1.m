clear all, close all, clc
% Schema explicite decentre
V = 0.1;     % vitesse d’advection
Tmax = 10;   % temps maximum
L = 5;       % longueur du domaine

N = 201;      % nb de points de discretisation en espace
dx = L/(N-1); % pas d’espace
s = 0.9;      % nombre de Courant
dt = s*dx/V;  % pas de temps

x = linspace(0,L,N)';
%condinit = @(x) (x>1.0) & (x<1.5);
condinit = @(x) sin(8*pi*x/L);

u = double(condinit(x)); % conversion booleen reel
tps = dt:dt:Tmax;
for t = tps
    uold = u;
    u(2:N) = uold(2:N)-s*(uold(2:N)-uold(1:N-1));
    u(1) = u(N);
end

uexact = condinit(x-V*t);
plot(x,u,'rx-',x,uexact,'b-'), grid on
title('Comparaison entre la solution exacte et approchee')
legend('Sol. approchee','Sol. exacte')
xlabel('x')
ylabel('Solution')