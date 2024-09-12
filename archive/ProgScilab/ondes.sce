////////////////////////////////////////////////////////////////
// Copyright G. Allaire, Juillet 2001, Janvier 2012
//
// Sur divers schemas de differences finies
// pour l'equation des ondes:
// comparaison des schemas
//
////////////////////////////////////////////////////////////////
//
// equation des ondes
// u,tt - a^2 u,xx = 0
//
////////////////////////////////////////////////////////////////
//
lg = 1. ; // lg = longueur du domaine
dx = 0.01 ; // dx = pas d'espace
nx = lg/dx ; // nx = nombre de mailles
a = 1. ; // a = vitesse
cfl = 0.9 ; // cfl
dt = dx*cfl/a ; // dt = pas de temps
tfinal = 5. ; // tfinal = temps de simulation
nt = int(tfinal/dt) ; // nt = nombre de pas de temps effectues
//
// initialisation: 
//
x=zeros(1,nx) ;
u0=zeros(1,nx) ;
u1=zeros(1,nx) ;
for i=1:nx
  x(i) = (i-1)*dx ;
  u0(i) = sin(x(i)*2*%pi) ;
  u1(i) = 0. ;
end
w=u0 ;
v=u0 ;
u=u0+dt*u1 ;
up=u0 ;
um=u0 ;
w4=u0 ;
v4=u0 ;
u4=u0+dt*u1 ;
up4=u0 ;
um4=u0 ;
uexacte=u0 ;

tics=[4,10,4,10];
plotframe([0.,-1.1,lg,1.1],tics);
plot2d(x,u0,1,"000")
xtitle ('donnee initiale' ,' ',' ') ;


halt() ;

// matrice du schema implicite
theta = 0.25 ;
mat=zeros(nx,nx) ;

for i=2:nx-1
  mat(i,i) = 1. + 2*theta*a*a*dt*dt/(dx*dx) ;
  mat(i,i+1) =  -theta*a*a*dt*dt/(dx*dx) ;
  mat(i,i-1) = -theta*a*a*dt*dt/(dx*dx) ;
end
mat(1,1) = 1. + 2*theta*a*a*dt*dt/(dx*dx) ;
mat(1,2) =  -theta*a*a*dt*dt/(dx*dx) ;
mat(nx,nx-1) = -theta*a*a*dt*dt/(dx*dx);
mat(nx,nx) = 1. + 2*theta*a*a*dt*dt/(dx*dx) ;
smat = sparse(mat) ;
spcho = chfact(smat) ; // factorisation de Cholesky

////////////////////////////////////////////////////////////////
// u = schema centre explicite
// (u= u^n, v=u^{n-1}, w=u^{n+1})
// u4 = theta-schema implicite
////////////////////////////////////////////////////////////////
//
for n=1:nt
//
up = shiftp('+1',u) ;
um = shiftp('-1',u) ;

w=2*u-v + a*a*dt*dt/(dx*dx)*(up-2*u+um) ;
v=u ;
u=w ;

up4 = shiftp('+1',u4) ;
um4 = shiftp('-1',u4) ;
vp4 = shiftp('+1',v4) ;
vm4 = shiftp('-1',v4) ;

b = 2*u4-v4 + (1-2*theta)*a*a*dt*dt/(dx*dx)*(up4-2*u4+um4) + theta*a*a*dt*dt/(dx*dx)*(vp4-2*v4+vm4) ;
w4 = chsolve(spcho,b) ; // resolution du systeme lineaire
v4=u4 ;
u4=w4 ;

if modulo(n,5) == 0
        for i=1:nx
            uexacte(i) = 0.5*( sin((x(i)-a*n*dt)*2*%pi) + sin((x(i)+a*n*dt)*2*%pi) ) ;
        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema centre explicite")
        plot2d(x,u4,[3,3],"100","theta-schema implicite, theta=0.25")
        xtitle('schemas pour l equation des ondes',' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        for i=1:nx
            uexacte(i) = 0.5*( sin((x(i)-a*nt*dt)*2*%pi) + sin((x(i)+a*nt*dt)*2*%pi) ) ;
        end
        clf()
        plotframe([0.,-1.1,lg,1.1],tics);
        plot2d(x,uexacte,[1,1],"100","solution exacte")
        plot2d(x,u,[2,2],"100","schema centre explicite")
        plot2d(x,u4,[3,3],"100","theta-schema implicite, theta=0.25")
        xtitle('schemas pour l equation des ondes, temps=5',' ',' ');


