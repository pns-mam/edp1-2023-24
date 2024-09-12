////////////////////////////////////////////////////////////////
// Copyright G. Allaire, Juillet 2001, Janvier 2012
//
// Sur divers schemas de differences finies
// pour l'equation de convection-diffusion:
// simulations 2-D
//
////////////////////////////////////////////////////////////////
//
// equation de convection-diffusion
// u,t + a u,x- u,xx = 0
//
////////////////////////////////////////////////////////////////
//
lgx = 3. ; // lgx = longueur du domaine
lgy = 1. ; // lgy = hauteur du domaine
dx = 0.05 ; // dx = pas d'espace
dy = 0.05 ; // dy = pas d'espace
nx = lgx/dx ; // nx = nombre de mailles en x
ny = lgy/dy ; // ny = nombre de mailles en y 
cfl = 0.4 ; // cfl
ax = 1. ; // ax = vitesse selon x
ay = 0. ; // ay = vitesse selon y
nu = 0.001 ; // nu = coefficient de diffusion
dt = (cfl/nu)/(1/(dx*dx)+1/(dy*dy)) ; // dt = pas de temps
dt = min( dt , dx*0.4/ax ) ;
tfinal = 0.5*lgx/ax ;
nt = floor(tfinal/dt) ; // nt = nombre de pas de temps effectues
//
// initialisation
//
x=zeros(1,nx+1) ;
y=zeros(1,ny+1) ;
u0=zeros(nx+1,ny+1) ;
for i=1:nx+1
  x(i) = (i-1)*dx  ;
  for j=1:ny+1
      y(j) = (j-1)*dy ;
      u0(i,j) = max(0.,1.-10*(x(i)-0.5)**2-10*(y(j)-0.5)**2) ;
  end
end
u=u0 ;
up=u0 ;
um=u0 ;
uexacte=u0 ;

//tics=[4,16,4,10];
//plotframe([0.,0.,lgx,lgy],tics);
//contour2d(x,y,u0,9,1:9,"000")
//xtitle ('donnee initiale' ,' ',' ') ;

leg='x@y@u' ;
flag=[0,1,4] ;
ebox=[0.,lgx,0.,lgy,0.,1.] ;
//plot3d(x,y,u0,35,45,leg) ;
plot3d(x,y,u0,35,45,leg,flag,ebox) ;
titre =  msprintf('donnee initiale') ; 
xtitle(titre,' ',' ');
halt() ;
////////////////////////////////////////////////////////////////
// schema explicite: cfl=0.4, vitesse=1.
////////////////////////////////////////////////////////////////
//
for n=1:nt
//
un = shift2d('n',u) ;
us = shift2d('s',u) ;
ue = shift2d('e',u) ;
uo = shift2d('o',u) ;

u=u - ax*dt/(dx)*(u-uo) - ay*dt/(2*dy)*(un-us) + nu*dt/(dx*dx)*(ue+uo-2*u) + nu*dt/(dy*dy)*(un+us-2*u) ;

if modulo(n,5) == 0
        clf()
//        flag=[0,0,4] ;
        plot3d(x,y,u,35,45,leg,flag,ebox) ;
        titre =  msprintf('Schema explicite, vitesse=%.2f, temps=%.2f, nombre de pas de temps=%i',ax,n*dt,n) ; 
        xtitle(titre,' ',' ');
        xpause(100000) ;
end
//
end

// comparaison au temps final
        titre =  msprintf('Schema explicite, vitesse=%.2f, temps final=%.2f, nombre de pas de temps=%i',ax,tfinal,nt) ; 
        clf()
        plot3d(x,y,u,35,45,leg,flag,ebox) ;
        xtitle(titre,' ',' ');

//[xx,yy,uu]=genfac3d(x,y,u) ;
//plot3d(xx,yy,uu,35,45,leg,flag) ;



