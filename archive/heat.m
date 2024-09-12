
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Equation de la chaleur
% u,t - u,xx = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation
%
lg = 10. ;   % lg = demie longueur du domaine
dx = 0.05 ;  % dx = pas d'espace
nx = lg/dx ; % nx = demi nombre de mailles
cfl = 0.4 ;  % cfl
dt = dx^2*cfl ; % dt = pas de temps
nt = 200 ;   % nt = nombre de pas de temps effectues

x = linspace(-lg,lg,2*nx+1);
u0 = max(0,1-x.^2);

u = u0 ;

plot(x,u0,'b-'),grid on
xlabel('x')
ylabel('Solution')
title('Donnee initiale') ;

%% schema explicite: cfl=0.4

for n = 1:nt
    u = u + dt/(dx*dx)*([u(2:end) 0]+[0 u(1:end-1)]-2*u) ;
    if mod(n,5) == 0
        plot(x,u,'r-',x,u0,'b-'), grid on
        legend('schema explicite','donnee initiale')
        title('schema explicite, cfl=0.4');
        pause(0.1)
    end
end
% comparaison au temps final
uexacte = zeros(1,2*nx+1) ;
noyau = zeros(1,2*nx+1) ;
for i=1:2*nx+1
    noyau(i) = exp(-((i-1)*dx)^2/(4*nt*dt)) ;
    if noyau(i) < 1.e-14
        kmax = i ;
        break
    end
end
for i=1:2*nx+1
    jmin = max(1, i-kmax) ;
    jmax = min(2*nx+1,i+kmax) ;
    for j=jmin:jmax
        k =1 + abs(i-j) ;
        uexacte(i) = uexacte(i) + u0(j)*dx*noyau(k) ;
    end
    uexacte(i) = uexacte(i)/sqrt(4*pi*nt*dt) ;
end

plot(x,u,'ro-',x,uexacte,'bx-'), grid on
legend('schema explicite','solution exacte')
title('schema explicite, cfl=0.4');

%% schema explicite: cfl=0.51
cfl = 0.501 ;
dt = dx*dx*cfl ;
u = u0 ;
nt = 1000 ;
for n = 1:nt
    u = u + dt/(dx*dx)*([u(2:end) 0]+[0 u(1:end-1)]-2*u) ;
    if mod(n,10) == 0
        plot(x,u,'r-',x,u0,'b-'), grid on
        legend('schema explicite','donnee initiale')
        title('schema explicite, cfl=0.51');
        pause(0.1)
    end
end

% comparaison au temps final
uexacte = zeros(1,2*nx+1) ;
noyau = zeros(1,2*nx+1) ;
for i=1:2*nx+1
    noyau(i) = exp(-((i-1)*dx)^2/(4*nt*dt)) ;
    if noyau(i) < 1.e-14
        kmax = i ;
        break
    end
end
for i=1:2*nx+1
    jmin = max(1, i-kmax) ;
    jmax = min(2*nx+1,i+kmax) ;
    for j=jmin:jmax
        k =1 + abs(i-j) ;
        uexacte(i) = uexacte(i) + u0(j)*dx*noyau(k) ;
    end
    uexacte(i) = uexacte(i)/sqrt(4*pi*nt*dt) ;
end

plot(x,u,'ro-',x,uexacte,'bx-'), grid on
legend('schema explicite','solution exacte')
title('schema explicite, cfl=0.51');
