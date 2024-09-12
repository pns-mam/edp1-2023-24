clear all, close all, clc
lg = 10. ;   % lg = demie longueur du domaine
dx = 0.05 ;  % dx = pas d'espace
nx = lg/dx ; % nx = demi nombre de mailles
cfl = 5 ;    % cfl
dt = dx^2*cfl ; % dt = pas de temps
nt = 150 ;   % nt = nombre de pas de temps effectues

x = linspace(-lg,lg,2*nx+1); %
u0 = max(0,1-x.^2);

A = diag((1+2*cfl)*ones(2*nx+1,1))+diag(-cfl*ones(2*nx,1),1)+...
    diag(-cfl*ones(2*nx,1),-1);
u = u0' ;

for n = 1:nt
    u = A\u ;
    if mod(n,5) == 0
        plot(x,u,'r-',x,u0,'b-'), grid on
        legend('schema implicite','donnee initiale')
        title('schema implicite, cfl=2');
        pause(0.1)
    end
end

% calcul de la solution exacte: convolution avec le noyau gaussien
uexacte = zeros(2*nx+1,1) ;
noyau = exp(-((0:2*nx)*dx).^2/(4*nt*dt)) ;
kmax = min(find(noyau<1.e-14));
for i=1:2*nx+1
    jmin = max(1, i-kmax) ;
    jmax = min(2*nx+1,i+kmax) ;
    uexacte(i) = sum(u0(jmin:jmax)*dx.*noyau(1 + abs(i-(jmin:jmax))))/sqrt(4*pi*nt*dt);
end

% comparaison au temps final
plot(x,u,'r-',x,uexacte,'b-'), grid on
legend('schema implicite','solution exacte')
title('schema implicite, cfl=2');
norm(u-uexacte,'inf')
