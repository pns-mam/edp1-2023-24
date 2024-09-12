////////////////////////////////////////////////////////////////
// Copyright G. Allaire, Juillet 2001, Janvier 2012
//
// Sur divers schemas de differences finies
// pour l'equation de la chaleur:
// mise en evidence d'instabilite, convergence
// trois schemas: centre instable, explicite, implicite
//
////////////////////////////////////////////////////////////////
//
// equation de la chaleur
// u,t - u,xx = 0
//
////////////////////////////////////////////////////////////////
//
lg = 10. ; // lg = demie longueur du domaine
dx = 0.05 ; // dx = pas d'espace
nx = lg/dx ; // nx = demi nombre de mailles
cfl = 0.4 ; // cfl
dt = dx*dx*cfl ; // dt = pas de temps
nt = 500 ; // nt = nombre de pas de temps effectues
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
// schema explicite: cfl=0.4
////////////////////////////////////////////////////////////////
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u=u + dt/(dx*dx)*(up+um-2*u) ;

if modulo(n,5) == 0
        clf()
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite, cfl=0.4',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        noyau = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            noyau(i) = exp(-((i-1)*dx)**2/(4*nt*dt)) ;
            if noyau(i) < 1.e-14 then
            kmax = i ;
            break
            end
        end
        for i=1:2*nx+1
            jmin = max(1, i-kmax) ;
            jmax = min(2*nx+1,i+kmax) ;
            for j=jmin:jmax
                 k =1 + int(abs(i-j)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau(k) ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf()
        plotframe([-lg,-0.1,lg,1.1],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema explicite, cfl=0.4, 500 pas de temps',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema explicite: cfl=0.51
////////////////////////////////////////////////////////////////
cfl = 0.51 ;
dt = dx*dx*cfl ;
u = u0 ;
nt=180 ;
//
for n=1:nt
//
up = shift('+1',u) ;
um = shift('-1',u) ;

u = u + dt/(dx*dx)*(up+um-2*u) ;

if modulo(n,5) == 0
        clf()
        plotframe([-lg,-0.8,lg,1.6],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema explicite, cfl=0.51',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        noyau = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            noyau(i) = exp(-((i-1)*dx)**2/(4*nt*dt)) ;
            if noyau(i) < 1.e-14 then
            kmax = i ;
            break
            end
        end
        for i=1:2*nx+1
            jmin = max(1, i-kmax) ;
            jmax = min(2*nx+1,i+kmax) ;
            for j=jmin:jmax
                 k =1 + int(abs(i-j)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau(k) ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf()
        plotframe([-lg,-0.8,lg,1.6],tics);
        plot2d(x,u,[1,1],"100","schema explicite")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema explicite, cfl=0.51, 180 pas de temps',' ',' ');

halt() ;
////////////////////////////////////////////////////////////////
// schema centre instable
////////////////////////////////////////////////////////////////
cfl = 0.1 ;
dt = dx*dx*cfl ;
u = u0 ; // solution au temps n+1
u1= u0 ; // solution au temps n
u2 = u0 ; // solution au temps n-1
nt=25 ;
//
// premier pas de temps
//
up = shift('+1',u0) ;
um = shift('-1',u0) ;

u1 = u0 + dt/(dx*dx)*(up+um-2*u0) ;

//
// pas de temps suivants
//
for n=2:nt
//
up = shift('+1',u1) ;
um = shift('-1',u1) ;

u = u2 + 2*dt/(dx*dx)*(up+um-2*u1) ;
u2 = u1 ;
u1 = u ;

if modulo(n,1) == 0
        clf()
        plotframe([-lg,-1.1,lg,2.1],tics);
        plot2d(x,u,[1,1],"100","schema centre")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema centre, cfl=0.1',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        noyau = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            noyau(i) = exp(-((i-1)*dx)**2/(4*nt*dt)) ;
            if noyau(i) < 1.e-14 then
            kmax = i ;
            break
            end
        end
        for i=1:2*nx+1
            jmin = max(1, i-kmax) ;
            jmax = min(2*nx+1,i+kmax) ;
            for j=jmin:jmax
                 k =1 + int(abs(i-j)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau(k) ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf()
        plotframe([-lg,-1.1,lg,2.1],tics);
        plot2d(x,u,[1,1],"100","schema centre")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema centre, cfl=0.1, 25 pas de temps',' ',' ');

halt() ;

////////////////////////////////////////////////////////////////
// schema implicite: cfl=2.
////////////////////////////////////////////////////////////////
cfl = 2. ;
dt = dx*dx*cfl ;
u = u0 ;
nt=200 ;
mat=zeros(2*nx+1,2*nx+1) ;

for i=2:2*nx
  mat(i,i) = 1. + 2*dt/(dx*dx) ;
  mat(i,i+1) =  -dt/(dx*dx) ;
  mat(i,i-1) = -dt/(dx*dx) ;
end
mat(1,1) = 1. + 2*dt/(dx*dx) ;
mat(1,2) =  -dt/(dx*dx) ;
mat(2*nx+1,2*nx) = -dt/(dx*dx);
mat(2*nx+1,2*nx+1) = 1. + 2*dt/(dx*dx) ;
smat = sparse(mat) ;
spcho = chfact(smat) ; // factorisation de Cholesky

//
for n=1:nt
//
u = chsolve(spcho,u) ; // resolution du systeme lineaire

if modulo(n,5) == 0
        clf()
        plotframe([-lg,-0.8,lg,1.6],tics);
        plot2d(x,u,[1,1],"100","schema implicite")
        plot2d(x,u0,[2,2],"100","donnee initiale")
        xtitle('schema implicite, cfl=2.',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        uexacte = zeros(1,2*nx+1) ;
        noyau = zeros(1,2*nx+1) ;
        for i=1:2*nx+1
            noyau(i) = exp(-((i-1)*dx)**2/(4*nt*dt)) ;
            if noyau(i) < 1.e-14 then
            kmax = i ;
            break
            end
        end
        for i=1:2*nx+1
            jmin = max(1, i-kmax) ;
            jmax = min(2*nx+1,i+kmax) ;
            for j=jmin:jmax
                 k =1 + int(abs(i-j)) ;
                 uexacte(i) = uexacte(i) + u0(j)*dx*noyau(k) ;
            end
            uexacte(i) = uexacte(i)/sqrt(4*%pi*nt*dt) ;
        end
        clf()
        plotframe([-lg,-0.8,lg,1.6],tics);
        plot2d(x,u,[1,1],"100","schema implicite")
        plot2d(x,uexacte,[2,2],"100","solution exacte")
        xtitle('schema implicite, cfl=2., 200 pas de temps',' ',' ');


