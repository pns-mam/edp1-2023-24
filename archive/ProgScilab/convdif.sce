////////////////////////////////////////////////////////////////
// Copyright G. Allaire, Juillet 2001
//
// Sur divers schemas de differences finies
// pour l'equation de convection-diffusion:
// mise en evidence d'instabilite, convergence
//
////////////////////////////////////////////////////////////////
//
// equation de convection-diffusion
// u,t + a u,x- u,xx = 0
//
////////////////////////////////////////////////////////////////
//
lg = 10. ; // lg = demie longueur du domaine
dx = 0.05 ; // dx = pas d'espace
nx = lg/dx ; // nx = demi nombre de mailles
cfl = 0.4 ; // cfl
dt = dx*dx*cfl ; // dt = pas de temps
nt = 300 ; // nt = nombre de pas de temps effectues
a = 1. ; // a = vitesse
//
// initialisation
//
x=zeros(1,2*nx+1) ;
u0=zeros(1,2*nx+1) ;
for i=1:2*nx+1
  x(i) = (i-nx-1)*dx  ;
  u0(i) = max(0.,1.-x(i)**2) ;
end
u=u0 ;
up=u0 ;
um=u0 ;
uexacte=u0 ;

tics=[4,10,4,10];
plotframe([-lg,-0.1,lg,1.1],tics);
plot2d(x,u0,1,"000")
xtitle ('donnee initiale' ,' ',' ') ;


halt() ;
////////////////////////////////////////////////////////////////
// schema explicite: cfl=0.4, vitesse=1.
////////////////////////////////////////////////////////////////
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u - a*dt/(2*dx)*(up-um) + dt/(dx*dx)*(up+um-2*u) ;

if modulo(n,5) == 0
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite, cfl=0.4, V=1',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            for j=1:2*nx+1
                 noyau =  exp(-((i-j)*dx-nt*dt*a)**2/(4*nt*dt)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema explicite, cfl=0.4, V=1, 300 pas de temps',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema explicite: cfl=0.4, vitesse=10.
////////////////////////////////////////////////////////////////
a=10. ;
u=u0 ;
nt=300;
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u - a*dt/(2*dx)*(up-um) + dt/(dx*dx)*(up+um-2*u) ;

if modulo(n,5) == 0
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite, cfl=0.4, V=10',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            for j=1:2*nx+1
                 noyau =  exp(-((i-j)*dx-nt*dt*a)**2/(4*nt*dt)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema explicite, cfl=0.4, V=10, 300 pas de temps',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema explicite: cfl=0.4, vitesse=100.
////////////////////////////////////////////////////////////////
nt=30 ;
a=100. ;
u=u0 ;
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u - a*dt/(2*dx)*(up-um) + dt/(dx*dx)*(up+um-2*u) ;

if modulo(n,5) == 0
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite, cfl=0.4, V=100',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            for j=1:2*nx+1
                 noyau =  exp(-((i-j)*dx-nt*dt*a)**2/(4*nt*dt)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema explicite, cfl=0.4, V=100, 30 pas de temps',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema explicite: cfl=0.1, vitesse=100.
////////////////////////////////////////////////////////////////
nt=120 ;
a=100. ;
u=u0 ;
cfl = 0.1 ; 
dt = dx*dx*cfl ; 
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u - a*dt/(2*dx)*(up-um) + dt/(dx*dx)*(up+um-2*u) ;

if modulo(n,5) == 0
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite, cfl=0.1, V=100',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            for j=1:2*nx+1
                 noyau =  exp(-((i-j)*dx-nt*dt*a)**2/(4*nt*dt)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema explicite, cfl=0.1, V=100, 120 pas de temps',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema explicite decentre: cfl=0.1, vitesse=100.
////////////////////////////////////////////////////////////////
nt=120 ;
a=100. ;
u=u0 ;
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u - a*dt/(dx)*(u-um) + dt/(dx*dx)*(up+um-2*u) ;

if modulo(n,5) == 0
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite decentre")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite decentre, cfl=0.1, V=100',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            for j=1:2*nx+1
                 noyau =  exp(-((i-j)*dx-nt*dt*a)**2/(4*nt*dt)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf;
        plotframe([-lg,-0.4,lg,1.2],tics);
        plot2d(x,u,[1,1],"100","schema explicite decentre")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema explicite decentre, cfl=0.1, V=100, 120 pas de temps',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// convection pure
// schema explicite centre: vitesse=100.
////////////////////////////////////////////////////////////////
nt=150 ;
a=100. ;
u=u0 ;
dt = 0.9*dx/a ;
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u - a*dt/(2*dx)*(up-um) ;

if modulo(n,5) == 0
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite centre")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite centre, V=100',' ',' ');
        xpause(100000) ;
end
//
end

halt() ;
////////////////////////////////////////////////////////////////
// convection pure
// schema explicite decentre: vitesse=100.
////////////////////////////////////////////////////////////////
nt=150 ;
a=100. ;
u=u0 ;
dt = 0.9*dx/a ;
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u - a*dt/(dx)*(u-um) ;

if modulo(n,5) == 0
        clf;
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite decentre")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite decentre, V=100',' ',' ');
        xpause(100000) ;
end
//
end

