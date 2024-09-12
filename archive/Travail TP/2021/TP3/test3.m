clear all, close all, clc
% Synthese: les deux schemas avec differentes conditions initiales
% On dessine la solution a chaque pas de temps
V = 0.1;     % vitesse d’advection
Tmax = 10;   % temps maximum
L = 5;       % longueur du domaine

N = 201;      % nb de points de discretisation en espace
dx = L/(N-1); % pas d’espace
s = 0.8;      % nombre de Courant
dt = s*dx/V;  % pas de temps
methode = 2;  % 1 = Explicite decentre, 2 = Lax-Wendroff
cinit = 1;    % 1 = creneau, 2 = sinus

x = linspace(0,L,N)';
if cinit == 1
    condinit = @(x) (x>1.0) & (x<1.5);
else
    condinit = @(x) sin(8*pi*x/L);
end
u = double(condinit(x));

tps = dt:dt:Tmax;
for t = tps
    uold = u;
    uexact = condinit(x-V*t);
    if methode == 1
        u(2:N) = uold(2:N)-s*(uold(2:N)-uold(1:N-1));
    else
        u(N) = uexact(N);
        u(2:N-1) = uold(2:N-1)-s*(uold(3:N)-uold(1:N-2))/2 +...
            s^2*(uold(3:N)-2*uold(2:N-1)+uold(1:N-2))/2 ;
    end
    u(1) = u(N);
    plot(x,u,'rx-',x,uexact,'bo-'), grid on
    title('Comparaison entre la solution exacte et approchee')
    legend('Sol. approchee','Sol. exacte')
    xlabel('x')
    ylabel('Solution')
    pause(0.1)
end